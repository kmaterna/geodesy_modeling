
from tectonic_utils.geodesy import haversine
import numpy as np
import os
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb.fault_slip_triangle import triangle_okada
from elastic_stresses_py import PyCoulomb
import elliptical_utilities  # local import
import inversion_utilities


default_params = PyCoulomb.configure_calc.get_lightweight_config_params(mu=30e9, lame1=30e9, B=0)


def elastic_model(param_vector, data_1d, cart_disp_points, faults, configs):
    """
    The forward model. The mesh geometry is known separately.  We take data points and a parameter vector,
    and convert them into more data points.
    Param vector must be 2n+3 parameters long -- n slip, n depth, 3 planar fit

    :param param_vector: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param data_1d: 1d_insar_object, which contains the look vector information for each pixel
    :param cart_disp_points: used only for location of points, in cartesian coordinates
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :param configs: dictionary of configuration parameters for the experiment
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    fault_list = []
    a, b, c = param_vector[-3], param_vector[-2], param_vector[-1]  # unpack parameters of the planar fit and refpix
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        surface_slip, bottom_depth = param_vector[j], param_vector[j+configs["num_faults"]]  # for real faults. KEY LINE
        surface_slip = surface_slip / 100  # converting from cm to m
        used_depths, widths, modeled_slip = elliptical_utilities.get_tuned_depth_arrays(configs['top_depth'],
                                                                                        configs['max_disloc_depth'],
                                                                                        configs["sampling_interval"],
                                                                                        surface_slip,
                                                                                        bottom_depth)
        for mslip, mwidth, depth in zip(modeled_slip, widths, used_depths):   # create all the source patches
            one_fault = fso.FaultSlipObject(strike=fault.strike, dip=89.99, length=fault.length,
                                            width=mwidth, lon=fault.lon,
                                            lat=fault.lat, depth=depth, rake=180, slip=mslip)
            # Change the slip to elliptical slip
            fault_list.append(one_fault)
    source_object = fso.fault_object_to_coulomb_fault(fault_list, zerolon_system=configs["zerolon"],
                                                      zerolat_system=configs["zerolat"])
    inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                    zerolon=configs['zerolon'],
                                                                                    zerolat=configs['zerolat'],
                                                                                    bbox=configs["bbox"])
    model_disp_points = triangle_okada.compute_cartesian_def_tris(inputs, default_params, cart_disp_points)  # run okada
    insar_1d_model = inversion_utilities.project_disp_points_into_insar1d(model_disp_points, data_1d)  # cart. space
    # Implement a planar fit to the data, which has three planar parameters for x and y and a constant offset
    insar_1d_model = insar_1d_model.subtract_ramp(a, b, c)
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm, in cartesian coordinates


def set_up_initial_params_and_bounds_full_vector(configs, arguments):
    """
    Here we construct an initial parameter vector of length 2N + 3.
    The last three are plane and refpix parameters.
    The structure of the parameter vector is [N*slip_m, N*depth_km, a, b, refpix_offset]
    The structure of lower bounds is [2N+3]
    The structure of upper bounds is [2N+3]
    """
    # Set up constraints: lower bound on slip = from fieldwork and creepmeters
    numfaults = configs["num_faults"]
    lower_bound = np.zeros((numfaults*2 + 3, ))  # lower bound on slip and depth starts at zero
    upper_bound = np.zeros((numfaults*2 + 3,))

    _, fieldwork_data = np.loadtxt(configs['fieldfile'], unpack=True)
    fieldwork_data = np.multiply(fieldwork_data, 0.1)  # convert mm to cm

    # Setting proper bounds
    lower_bound[0:numfaults] = fieldwork_data  # lower bound on slip from field data
    upper_bound[0:numfaults] = 6  # upper bound on slip is 6 cm
    lower_bound[numfaults:2*numfaults] = 0  # lower bound on depth is zero
    upper_bound[numfaults:2*numfaults] = 5  # upper bound on depth is 5 km
    lower_bound[-3], lower_bound[-2], lower_bound[-1] = -50, -50, -10.0  # parameters for plane and reference pixel
    upper_bound[-3], upper_bound[-2], upper_bound[-1] = 50, 50, 10.0

    if configs['starting_point'] is not None:
        param0 = np.loadtxt(configs['starting_point'])
    else:
        param0 = []   # param vector = [slip slip slip ..... slip depth depth depth..... a, b, c]
        for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
            slip_guess = np.max([2.0, lower_bound[i]])  # the initial guess is either 20 mm or the minimum of range
            param0.append(slip_guess)
        for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
            param0.append(2.0)
        param0 = param0 + [0, 0, 0]  # plane, plane, and reference pixel

    # Set up the original parameter vector and the upper bounds and lower bounds
    x_scale = []  #
    for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
        x_scale.append(1.0)  # scale is mm
    for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
        x_scale.append(1.0)
    x_scale = x_scale + [5.0, 5.0, 1.0]
    value_and_bounds = np.vstack((param0, lower_bound, upper_bound)).T
    print("Writing initial guess in %s" % os.path.join(arguments.output, 'initial_guess.txt'))
    np.savetxt(os.path.join(arguments.output, 'initial_guess.txt'), value_and_bounds)
    return param0, lower_bound, upper_bound, x_scale


def set_up_initial_params_and_bounds_const_stress(configs, arguments):
    """
    Here we construct an initial parameter vector of length N + 4.
    N slip values, one strain drop.
    The last three are plane and refpix parameters.
    The structure of the parameter vector is [N*slip_m, strain, a, b, refpix_offset]
    The structure of lower bounds is [N+4]
    The structure of upper bounds is [N+4]
    """
    # Set up constraints: lower bound on slip = from fieldwork and creepmeters
    numfaults = configs["num_faults"]
    lower_bound = np.zeros((numfaults + 4, ))  # lower bound on slip starts at zero
    upper_bound = np.zeros((numfaults + 4,))

    _, fieldwork_data = np.loadtxt(configs['fieldfile'], unpack=True)
    fieldwork_data = np.multiply(fieldwork_data, 0.1)  # convert mm to cm

    # Setting proper bounds
    lower_bound[0:numfaults] = fieldwork_data  # lower bound on slip from field data
    upper_bound[0:numfaults] = 6  # upper bound on slip is 6 cm
    lower_bound[-4] = 0.02   # strain (e-5)
    upper_bound[-4] = 10  # strain (e-5)
    lower_bound[-3], lower_bound[-2], lower_bound[-1] = -50, -50, -10.0  # parameters for plane and reference pixel
    upper_bound[-3], upper_bound[-2], upper_bound[-1] = 50, 50, 10.0

    if configs['starting_point'] is not None:
        param0 = np.loadtxt(configs['starting_point'])
    else:
        param0 = []   # param vector = [slip slip slip ..... slip, strain, a, b, c]
        for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
            slip_guess = np.max([2.0, lower_bound[i]])  # the initial guess is either 20 mm or the minimum of range
            param0.append(slip_guess)
        param0 = param0 + [0, 0, 0, 0]  # strain, plane, plane, and reference pixel

    # Set up the original parameter vector and the upper bounds and lower bounds
    x_scale = []  #
    for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
        x_scale.append(1.0)  # scale is mm
    x_scale = x_scale + [1.0, 5.0, 5.0, 1.0]
    value_and_bounds = np.vstack((param0, lower_bound, upper_bound)).T
    print("Writing initial guess in %s" % os.path.join(arguments.output, 'initial_guess.txt'))
    np.savetxt(os.path.join(arguments.output, 'initial_guess.txt'), value_and_bounds)
    return param0, lower_bound, upper_bound, x_scale


def modify_cov_to_split_across_fault(matrix, data, faults):
    """
    Create a "lobotomized" version of the classic insar covariance matrix.
    Pixel-pairs close to the fault trace but on opposite sides of it can be fit more with a standard lowest RMS
    misfit rather than respecting the spatial structure of atmospheric noise.
    We want these pixels to be treated as signal.

    :param matrix: covariance matrix determined from noise data away from the fault
    :param data: insar1d object that has the lon/lat of pixel values
    :param faults: list of fault_slip_objects, used for fault traces
    :return: another covariance matrix
    """
    # Define parameters
    fault_strike = 315  # approximate fault strike
    jitter = 30  # adding jitter * identity matrix to keep covariance matrix stable. Results very sensitive to jitter
    close_distance = 1.5  # if pixels are closer than this distance, the covariance between them will be reduced
    lobotomy_value = 1  # the lobotomized pixels will have the following small covariance, in mm.
    # This should be smaller than the given covariances to reduce the weighting on these off-diagonals.

    fault_lons, fault_lats = [], []
    for i in faults:
        lons, lats = i.get_updip_corners_lon_lat()
        fault_lons.append(lons[0])
        fault_lons.append(lons[1])
        fault_lats.append(lats[0])
        fault_lats.append(lats[1])

    is_NE = []
    distances = []

    # Create a vector of whether each point is on the NE or SW side of the fault
    for i in range(len(data.LOS)):
        pt_lon, pt_lat = data.lon[i], data.lat[i]
        smallest_distance = 100  # arbitrary large number of km, will be reduced over iterations
        proper_heading = 0
        for x, y in zip(fault_lons, fault_lats):
            small_distance = haversine.distance((y, x), (pt_lat, pt_lon))
            heading = haversine.calculate_initial_compass_bearing((y, x), (pt_lat, pt_lon))
            if small_distance < smallest_distance:
                smallest_distance = small_distance
                proper_heading = heading
        distances.append(smallest_distance)  # in km
        if 0 < proper_heading < fault_strike-180 or fault_strike < proper_heading < 360:
            is_NE.append(1)  # if the point is northeast of the fault
        else:
            is_NE.append(0)  # if the point is southwest of the fault

    # 4) Eigen sanity (Hermitian eigvals are cheap & reliable)
    w = np.linalg.eigvalsh(matrix)
    print("eig min/max =", w[0], w[-1])
    # 5) Condition number proxy
    cond = w[-1] / max(w[0], 1e-300)
    print("Before changes: cond(eig) ≈", cond)

    m, n = np.shape(matrix)
    counter = 0
    for i in range(m):
        for j in range(i):
            pt1 = (data.lat[i], data.lon[i])
            pt2 = (data.lat[j], data.lon[j])
            distance = haversine.distance(pt1, pt2)
            if is_NE[i] + is_NE[j] == 1:  # if the two pixels are on opposite sides of the fault
                if distance < close_distance:  # if they are really close together, 1.5 km apart
                    matrix[i][j] = matrix[j][i] = lobotomy_value  # change the values to a small number
                    counter = counter + 1
    print(f"{counter} pixels out of {np.size(matrix)} pixels changed")

    # Ensuring the matrix is stable and symmetric
    matrix = matrix + jitter * np.eye(matrix.shape[0])
    matrix = 0.5 * (matrix + matrix.T)  # line is not exactly necessary, but it ensures matrix is symmetric

    # 4) Eigen sanity (Hermitian eigvals are cheap & reliable)
    w = np.linalg.eigvalsh(matrix)
    print("eig min/max =", w[0], w[-1])
    # 5) Condition number proxy
    cond = w[-1] / max(w[0], 1e-300)
    print("After mods and jitter: cond(eig) ≈", cond)

    return matrix
