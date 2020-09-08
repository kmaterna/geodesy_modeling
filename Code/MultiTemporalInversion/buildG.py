import numpy as np
import matplotlib.pyplot as plt
import sys
import slippy.xyz2geo as plotting_library
import slippy.basis
import slippy.patch
import slippy.tikhonov
import slippy.gbuild
import scipy.optimize
import scipy.linalg
import slippy.io


def reg_nnls(Gext, dext):
    return scipy.optimize.nnls(Gext, dext)[0]


def G_with_smoothing(G, L, alpha, d):
    # Add smoothing and other regularization
    # Expand the data to match.
    dext = np.concatenate((d, np.zeros(L.shape[0])))
    Gext = np.vstack((G, L))
    if alpha > 0:  # Minimum norm solution. Aster and Thurber, Equation 4.5.
        alphaI = alpha * np.identity(np.shape(Gext)[1])
        alphaI[-1, -1] = 0;  # for leveling, we don't want the offset term to be constrained.
        zero_vector = np.zeros((np.shape(Gext)[1],))
        Gext = np.vstack((Gext, alphaI))
        dext = np.concatenate((dext, zero_vector))
    return Gext, dext;


def input_gps_file(filename):
    obs_disp_f = np.zeros((0,))
    obs_sigma_f = np.zeros((0,))
    obs_pos_geo_f = np.zeros((0, 3))
    obs_basis_f = np.zeros((0, 3))

    gps_input = slippy.io.read_gps_data(filename);
    Ngps = len(gps_input[0])
    obs_gps_pos_geo = gps_input[0]
    obs_gps_disp = gps_input[1]
    obs_gps_sigma = gps_input[2]
    obs_gps_basis = slippy.basis.cardinal_basis((Ngps, 3))

    obs_disp_fi = obs_gps_disp.reshape((Ngps * 3,))
    obs_sigma_fi = obs_gps_sigma.reshape((Ngps * 3,))
    obs_basis_fi = obs_gps_basis.reshape((Ngps * 3, 3))
    obs_pos_geo_fi = obs_gps_pos_geo[:, None, :].repeat(3, axis=1).reshape((Ngps * 3, 3))

    obs_disp_f = np.concatenate((obs_disp_f, obs_disp_fi), axis=0)
    obs_sigma_f = np.concatenate((obs_sigma_f, obs_sigma_fi), axis=0)
    obs_basis_f = np.concatenate((obs_basis_f, obs_basis_fi), axis=0)
    obs_pos_geo_f = np.concatenate((obs_pos_geo_f, obs_pos_geo_fi), axis=0)

    print("Reading %d gps data from %s " % (len(gps_input[0]), filename));

    return [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ngps];


def input_insar_file(filename):
    obs_disp_f = np.zeros((0,))
    obs_sigma_f = np.zeros((0,))
    obs_pos_geo_f = np.zeros((0, 3))
    obs_basis_f = np.zeros((0, 3))

    insar_input = slippy.io.read_insar_data(filename)
    Ninsar = len(insar_input[0])
    obs_insar_pos_geo = insar_input[0]
    obs_insar_disp = insar_input[1]
    obs_insar_sigma = insar_input[2]
    obs_insar_basis = insar_input[3]

    obs_disp_f = np.concatenate((obs_disp_f, obs_insar_disp), axis=0)
    obs_sigma_f = np.concatenate((obs_sigma_f, obs_insar_sigma), axis=0)
    obs_basis_f = np.concatenate((obs_basis_f, obs_insar_basis), axis=0)
    obs_pos_geo_f = np.concatenate((obs_pos_geo_f, obs_insar_pos_geo), axis=0)

    print("Reading %d insar data from %s " % (len(insar_input[0]), filename));

    return [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ninsar];


def beginning_calc(config):
    # Read Faults into a list of faults (each fault being a dict)
    fault_list = [];
    for key in config["faults"].keys():
        basis = []
        for b in ['basis1', 'basis2', 'basis3']:
            basis_i = config["faults"][key].pop(b, None)
            if basis_i is None:
                continue
            basis += [basis_i]
        fault_segment = {
            "strike": config["faults"][key]["strike"],
            "dip": config["faults"][key]["dip"],
            "length": config["faults"][key]["length"],
            "width": config["faults"][key]["width"],
            "seg_pos_geo": config["faults"][key]["position"],
            "Nlength": config["faults"][key]["Nlength"],
            "Nwidth": config["faults"][key]["Nwidth"],
            "penalty": config["faults"][key]["penalty"],
            "slip_basis": basis};
        fault_list.append(fault_segment)
        alpha = config['alpha']

    # Setting up the basemap before we begin (using the GPS as information)
    for dataset in config["data_files"].keys():
        if config["data_files"][dataset]["type"] == "gps":
            gps_input = slippy.io.read_gps_data(config["data_files"][dataset]["data_file"]);
            Ngps = len(gps_input[0])
            obs_gps_pos_geo = gps_input[0]
            obs_pos_geo_basemap = obs_gps_pos_geo[:, None, :].repeat(3, axis=1).reshape((Ngps * 3, 3))
            break;
    bm = plotting_library.create_default_basemap(obs_pos_geo_basemap[:, 0], obs_pos_geo_basemap[:, 1])
    number_of_datasets = len(config["data_files"].keys());

    # # ###################################################################
    # ### discretize the fault segments
    # ### create slip basis vectors for each patch
    # ### build regularization matrix
    # ###################################################################
    patches = [];  # a growing list of fault patches.
    slip_basis_f = np.zeros((0, 3))  # a growing list of basis functions for slip patches
    total_fault_slip_basis = [];
    patches_f = [];
    L_array = [];
    Ns_total = 0;

    # # Fault processing
    for fault in fault_list:
        # Convert fault to cartesian coordinates
        fault["seg_pos_cart"] = plotting_library.geodetic_to_cartesian(fault["seg_pos_geo"], bm)

        # Discretize fault segment
        seg = slippy.patch.Patch(fault["seg_pos_cart"],
                                 fault["length"], fault["width"],
                                 fault["strike"], fault["dip"])
        single_fault_patches = np.array(seg.discretize(fault["Nlength"], fault["Nwidth"]))
        Ns = len(single_fault_patches);

        # Create slip basis vectors
        Ds = len(fault["slip_basis"])  # the number of basis vectors for this slip patch.
        single_fault_slip_basis = np.array([fault["slip_basis"] for j in range(Ns)])
        if total_fault_slip_basis == []:
            total_fault_slip_basis = single_fault_slip_basis;
        else:
            total_fault_slip_basis = np.concatenate((total_fault_slip_basis, single_fault_slip_basis), axis=0)
        single_fault_silp_basis_f = single_fault_slip_basis.reshape((Ns * Ds, 3))

        # Packaging of slip_basis_f
        single_fault_patches_f = single_fault_patches[:, None].repeat(Ds, axis=1).reshape((Ns * Ds,))
        patches = np.concatenate((patches, single_fault_patches), axis=0)
        patches_f = np.concatenate((patches_f, single_fault_patches_f), axis=0)
        slip_basis_f = np.concatenate((slip_basis_f, single_fault_silp_basis_f), axis=0);

        ### build regularization matrix
        L = np.zeros((0, Ns * Ds))
        indices = np.arange(Ns * Ds).reshape((Ns, Ds))
        for i in range(Ds):
            connectivity = indices[:, i].reshape((fault["Nlength"], fault["Nwidth"]))
            Li = slippy.tikhonov.tikhonov_matrix(connectivity, 2, column_no=Ns * Ds)
            L = np.vstack((Li, L))

        L *= fault["penalty"]
        L_array.append(L)
        Ns_total = Ns_total + Ns

    L = scipy.linalg.block_diag(*L_array)  # Make the block diagonal matrix for tikhonov regularization

    # PARSE HOW MANY EPOCHS WE ARE USING
    # This will tell us how many epochs, model parameters total, and output files we need.
    n_epochs = 0;
    total_spans = [];
    span_output_files = [];
    for epoch in config["epochs"].keys():
        n_epochs = n_epochs + 1;
        total_spans.append(config["epochs"][epoch]["name"]);
        span_output_files.append(config["epochs"][epoch]["slip_output_file"]);
    n_model_params = np.shape(L)[0];
    print("Finding fault model for %d epochs " % n_epochs);
    print("Number of fault model parameters per epoch: %d" % n_model_params);
    print("Total number of model parameters: %d" % (n_model_params * n_epochs + number_of_datasets));

    # INITIALIZING THE LARGE MATRICIES
    d_total = np.zeros((0,));
    sig_total = np.zeros((0,));
    weight_total = np.zeros((0,));
    G_total = np.zeros((0, n_epochs * n_model_params));  # does not contain leveling offsets
    G_nosmooth_total = np.zeros((0, n_epochs * n_model_params));  # for computing predicted displacements
    nums_obs = [];
    pos_obs_list = [];
    pos_basis_list = [];

    # INITIAL DATA SCOPING: HOW MANY FILES WILL NEED TO BE READ?
    input_file_list = [];
    spans_list = [];
    data_type_list = [];
    strengths_list = [];
    signs_list = [];
    output_file_list = [];
    row_span_list = [];
    print("Available data indicated in json file: ")

    for data_file in config["data_files"].keys():
        this_data_spans = config["data_files"][data_file]["span"];  # span expected of all data files
        this_strength = config["data_files"][data_file]["strength"];  # strength expected of all data files
        input_file_list.append(config["data_files"][data_file]["data_file"]);  # infile expected of all data files
        output_file_list.append(config["data_files"][data_file]["outfile"]);  # outfile expected of all data files
        data_type_list.append(config["data_files"][data_file]["type"]);  # type expected of all data files
        spans_list.append(this_data_spans);
        strengths_list.append(this_strength);
        print(config["data_files"][data_file]["data_file"], " spans ", this_data_spans);
        if config["data_files"][data_file]["type"] == "leveling":
            signs_list.append(config["data_files"][data_file]["sign"]);
        else:
            signs_list.append(0.0);

    # START THE MAJOR INVERSION LOOP
    for filenum in range(len(input_file_list)):
        leveling_sign = signs_list[filenum];

        # INPUT STAGE
        filename = input_file_list[filenum];
        if data_type_list[filenum] == 'gps':
            [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ngps] = input_gps_file(filename);
            obs_weighting_f = (1 / strengths_list[filenum]) * np.ones((Ngps * 3,));
        elif data_type_list[filenum] == 'insar':
            [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ninsar] = input_insar_file(filename);
            obs_weighting_f = (1 / strengths_list[filenum]) * np.ones((Ninsar,));
        elif data_type_list[filenum] == "leveling":
            [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ninsar] = input_insar_file(filename);
            obs_weighting_f = (1 / strengths_list[filenum]) * np.ones((Ninsar,));
        else:
            print("ERROR! Unrecognized data type %s " % data_type_list[filenum]);
            continue;

        # Metadata collection about this input source
        pos_obs_list.append(obs_pos_geo_f);  # for writing the outputs later
        nums_obs.append(len(obs_disp_f));
        pos_basis_list.append(obs_basis_f);
        this_data_spans = spans_list[filenum];
        print("This data spans:", this_data_spans);

        # # Geodetic to Cartesian coordinates
        obs_pos_cart_f = plotting_library.geodetic_to_cartesian(obs_pos_geo_f, bm)

        # BUILD SYSTEM MATRIX
        ###################################################################
        G = slippy.gbuild.build_system_matrix(obs_pos_cart_f,
                                              patches_f,
                                              obs_basis_f,
                                              slip_basis_f,
                                              leveling=False)

        # ### weigh system matrix and data by the uncertainty
        # ###################################################################
        G /= obs_weighting_f[:, None]
        obs_disp_f /= obs_weighting_f
        G /= obs_sigma_f[:, None]
        obs_disp_f /= obs_sigma_f
        print("G:", np.shape(G));

        # regularization matrix (tikhonov regularization)
        print("L:", np.shape(L))

        # BUILD THE G MATRIX FOR THIS SET OF OBSERVATIONS
        G_ext, d_ext = G_with_smoothing(G, L, alpha, obs_disp_f);
        print("Gext (G,L,alpha):", np.shape(G_ext));

        # APPEND TO THE BIGGER MATRIX AND DATA VECTOR (WITH THE CORRECT SPANS)
        n_rows = len(d_ext);
        n_cols = n_model_params * n_epochs;
        G_rowblock_obs = np.zeros((n_rows, n_cols));  # one per observation
        G_nosmoothing_obs = np.zeros((np.shape(G)[0], n_cols));  # one per observation

        # Which spans does the data cover? This is like the "1" in SBAS
        # One new "rowblock" for each data file
        # Also build the G matrix that only has data (no smoothing)
        count = 0;
        for epoch in total_spans:
            col_limits = [count * n_model_params, (count + 1) * n_model_params];
            if epoch in this_data_spans:
                G_rowblock_obs[:, col_limits[0]:col_limits[1]] = G_ext;
                G_nosmoothing_obs[:, col_limits[0]:col_limits[1]] = G;
            count = count + 1;

        # APPEND TO G_TOTAL AND G_NOSMOOTH_TOTAL
        d_total = np.concatenate((d_total, d_ext));  # Building up the total data vector
        sig_total = np.concatenate((sig_total, obs_sigma_f));  # building up the total sigma vector
        weight_total = np.concatenate((weight_total, obs_weighting_f));  # building up the total weighting vector
        top_row_data = len(G_total);  # Where in the matrix is this data?
        bottom_row_data = len(G_total) + len(G);  # Where in the matrix is this data?
        row_span_list.append([top_row_data, bottom_row_data]);
        G_total = np.concatenate((G_total, G_rowblock_obs));
        print("  Adding %d lines " % len(G_rowblock_obs))
        G_nosmooth_total = np.concatenate((G_nosmooth_total, G_nosmoothing_obs));

    # ADDING THE COLUMNS FOR LEVELING OFFSETS TO G_TOTAL MATRIX (Corresponding to the data lines only)
    print("------\nBefore adding lines for leveling offsets, shape(G): ", np.shape(G_total));
    numrows = np.shape(G_total)[0];
    num_rows_nosmooth = np.shape(G_nosmooth_total)[0];
    count = 0;
    for filenum in range(len(input_file_list)):
        leveling_sign = signs_list[filenum];
        newcol = np.zeros((numrows, 1));
        for i in range(len(newcol)):
            if i >= row_span_list[filenum][0] and i < row_span_list[filenum][1]:
                newcol[i] = leveling_sign;
        G_total = np.hstack((G_total, newcol));
        print("Adding column for %s " % input_file_list[filenum])

        # ADDING THE COLUMNS TO G_NOSMOOTH_TOTAL
        newcol_nosmooth = np.zeros((num_rows_nosmooth, 1));
        for i in range(len(newcol_nosmooth)):
            if i >= count and i < count + (row_span_list[filenum][1] - row_span_list[filenum][0]):
                newcol_nosmooth[i] = leveling_sign;
        count = count + (row_span_list[filenum][1] - row_span_list[filenum][0]);
        G_nosmooth_total = np.hstack((G_nosmooth_total, newcol_nosmooth));

    plt.figure(figsize=(12, 8), dpi=300);
    plt.imshow(G_total, vmin=-0.02, vmax=0.02, aspect=1/10)
    plt.savefig("image_of_G.png");

    # INVERT BIG-G
    #   ### estimate slip and compute predicted displacement
    #   #####################################################################
    slip_f = reg_nnls(G_total, d_total)  # the model
    print("G_total:", np.shape(G_total));
    print("slip_f:", np.shape(slip_f))
    print("G_nosmooth_total:", np.shape(G_nosmooth_total));
    print(np.shape(G_nosmooth_total.dot(slip_f)))
    print(np.shape(sig_total))
    print(np.shape(weight_total));
    pred_disp_f = G_nosmooth_total.dot(slip_f) * sig_total * weight_total;

    total_cardinal_slip, leveling_offsets = parse_slip_outputs(slip_f, Ns_total, Ds, n_epochs, total_fault_slip_basis,
                                                               number_of_datasets);
    disp_models = parse_disp_outputs(pred_disp_f, nums_obs);

    # Defensive programming (Reporting errors)
    for i in range(len(leveling_offsets)):
        if signs_list[i] != 0:
            print("Leveling Offset for %s = %f m " % (input_file_list[i], leveling_offsets[i]));
            if abs(leveling_offsets[i]) < 0.0000001:
                print("WARNING: Leveling offset for %s close to zero. Consider a negative offset in G." % (input_file_list[i]));

    ### get slip patch data for outputs
    #####################################################################
    patches_pos_cart = [i.patch_to_user([0.5, 1.0, 0.0]) for i in patches]
    patches_pos_geo = plotting_library.cartesian_to_geodetic(patches_pos_cart, bm)
    patches_strike = [i.strike for i in patches]
    patches_dip = [i.dip for i in patches]
    patches_length = [i.length for i in patches]
    patches_width = [i.width for i in patches]

    # OUTPUT EACH SLIP INTERVAL
    for i in range(n_epochs):
        one_span = total_spans[i];
        slip_output_file = span_output_files[i];

        # ### write output
        # #####################################################################
        slippy.io.write_slip_data(patches_pos_geo,
                                  patches_strike, patches_dip,
                                  patches_length, patches_width,
                                  total_cardinal_slip[i], slip_output_file)
        print("Writing file %s " % slip_output_file);

    # OUTPUT EACH PREDICTED DISPLACEMENT FIELD
    for filenum in range(len(input_file_list)):

        if data_type_list[filenum] == 'gps':  # write GPS
            Ngps = int(nums_obs[filenum] / 3);
            pred_disp_gps = disp_models[filenum];
            pred_disp_gps = pred_disp_gps.reshape((Ngps, 3))
            slippy.io.write_gps_data(pos_obs_list[filenum][::3],
                                     pred_disp_gps, 0.0 * pred_disp_gps,
                                     output_file_list[filenum]);
            print("Writing file %s " % output_file_list[filenum]);

        elif data_type_list[filenum] == 'insar':
            Npts = nums_obs[filenum];
            pred_disp_insar = disp_models[filenum];
            slippy.io.write_insar_data(pos_obs_list[filenum],
                                       pred_disp_insar, 0.0 * pred_disp_insar,
                                       pos_basis_list[filenum],
                                       output_file_list[filenum])
            print("Writing file %s " % output_file_list[filenum]);

        elif data_type_list[filenum] == 'leveling':
            Npts = nums_obs[filenum];
            pred_disp_leveling = disp_models[filenum];
            slippy.io.write_insar_data(pos_obs_list[filenum],
                                       pred_disp_leveling, 0.0 * pred_disp_leveling,
                                       pos_basis_list[filenum],
                                       output_file_list[filenum])
            print("Writing file %s " % output_file_list[filenum]);

    return;


def parse_slip_outputs(slip_f, Ns_total, Ds, n_epochs, total_fault_slip_basis, number_of_datasets):
    print("Parsing slip outputs from %d epochs " % (n_epochs));
    count = 0;
    total_cardinal_slip = [];
    for i in range(n_epochs):
        n_params = int((len(slip_f) - number_of_datasets) / n_epochs);  # num fault params; ex: 200
        start = count;
        finish = count + n_params;
        slip_during_epoch = slip_f[start:finish];
        count = count + n_params;

        slip = slip_during_epoch.reshape(
            (Ns_total, Ds))  # THIS ASSUMES ALL FAULTS HAVE THE SAME NUMBER OF BASIS VECTORS
        cardinal_slip = slippy.basis.cardinal_components(slip, total_fault_slip_basis)
        total_cardinal_slip.append(cardinal_slip);

    leveling_offsets = slip_f[-number_of_datasets:];

    return total_cardinal_slip, leveling_offsets;


def parse_disp_outputs(pred_disp_f, num_obs):
    count = 0;
    disp_segments = [];
    for i in range(len(num_obs)):
        disp_segment = pred_disp_f[count:count + num_obs[i]]
        disp_segments.append(disp_segment);
        count = count + num_obs[i];
    return disp_segments;
