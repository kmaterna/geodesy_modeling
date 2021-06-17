import numpy as np
import matplotlib.pyplot as plt
import json
import slippy.xyz2geo as plotting_library
import slippy.basis
import slippy.patch
import slippy.tikhonov
import slippy.gbuild
import scipy.optimize
import scipy.linalg
import slippy.io
from . import resolution_tests


def reg_nnls(Gext, dext):
    return scipy.optimize.nnls(Gext, dext)[0]


def G_with_smoothing(G, L, alpha, d):
    """
    L: Add smoothing regularization (tikhonov regularization)
    Alpha: Add minimum-norm regularization (Aster and Thurber, Equation 4.5)
    d: Expand the data vector to match the new size of G
    """
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
    """
    obs_disp_f: list of obs (8 stations --> 24 obs)
    obs_sigma_f: list of sigmas (8 stations --> 24 sigmas)
    obs_basis_f: 3-vector basis component (e, n, or u) for each component of each obs (8 stations --> 24 vectors)
    obs_pos_geo_f: 3-vector llh for each component of each obs (8 stations --> 24 llh)
    Ngps: int (8 for 8 stations)
    """
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


def input_faults(config):
    # Read Faults into a list of faults (each fault being a dict), starting from a filename for each fault
    print("Parsing faults from config file");
    fault_list = [];
    for i, key in enumerate(config["faults"].keys()):

        config_file = open(config["faults"][key]['filename'], 'r');
        fault_i = json.load(config_file);
        fault_i = fault_i[key]
        basis = []
        for b in ['basis1', 'basis2', 'basis3']:
            basis_i = fault_i.pop(b, None)
            if basis_i is None:
                continue
            basis += [basis_i]
        fault_segment = {
            "strike": fault_i["strike"],
            "dip": fault_i["dip"],
            "length": fault_i["length"],
            "width": fault_i["width"],
            "seg_pos_geo": fault_i["position"],
            "Nlength": fault_i["Nlength"],
            "Nwidth": fault_i["Nwidth"],
            "penalty": config["faults"][key]["penalty"],
            "slip_basis": basis,
            "name": i
        };
        fault_list.append(fault_segment)

    return fault_list;


def graph_big_G(config, G):
    # Show big-G matrix for all times, all data
    plt.figure(figsize=(12, 8), dpi=300);
    plt.imshow(G, vmin=-0.02, vmax=0.02, aspect=1/10);
    plt.savefig(config['output_dir']+"/image_of_G.png");
    return;


def parse_slip_outputs(slip_f, Ns_total, Ds, n_epochs, total_fault_slip_basis, number_of_datasets):
    """ simple utility function to reduce lines in the main program"""
    print("Parsing slip outputs from %d epochs " % n_epochs);
    count = 0;
    total_cardinal_slip = [];
    for i in range(n_epochs):
        n_params = int((len(slip_f) - number_of_datasets) / n_epochs);  # num fault params; ex: 200
        start = count;
        finish = count + n_params;
        slip_during_epoch = slip_f[start:finish];
        count = count + n_params;

        slip = slip_during_epoch.reshape((Ns_total, Ds))  # ASSUMES ALL FAULTS HAVE SAME NUMBER OF BASIS VECTORS
        cardinal_slip = slippy.basis.cardinal_components(slip, total_fault_slip_basis)
        total_cardinal_slip.append(cardinal_slip);

    leveling_offsets = slip_f[-number_of_datasets:];

    return total_cardinal_slip, leveling_offsets;


def parse_disp_outputs(pred_disp_f, num_obs):
    """ simple utility function """
    count = 0;
    disp_segments = [];
    for i in range(len(num_obs)):
        disp_segment = pred_disp_f[count:count + num_obs[i]]
        disp_segments.append(disp_segment);
        count = count + num_obs[i];
    return disp_segments;


def beginning_calc(config):

    with open(config['output_dir']+'/config.json', 'w') as fp:
        json.dump(config, fp, indent="  ");   # save copy of config file in outdir, for record-keeping

    fault_list = input_faults(config);
    alpha = config['alpha']  # a parameter to produce Minimum norm solution (optional)

    # Setting up the basemap before we begin (using the first dataset as information)
    number_of_datasets = len(config["data_files"].keys());   # integer
    first_dataset = list(config["data_files"].keys())[0]
    if config["data_files"][first_dataset]["type"] == "gps":
        first_input = slippy.io.read_gps_data(config["data_files"][first_dataset]["data_file"]);  # gps
    else:
        first_input = slippy.io.read_insar_data(config["data_files"][first_dataset]["data_file"]);   # lev or insar
    obs_pos_geo = first_input[0]
    obs_pos_geo_basemap = obs_pos_geo[:, None, :].repeat(3, axis=1).reshape((len(first_input[0]) * 3, 3))  # reshape llh
    bm = plotting_library.create_default_basemap(obs_pos_geo_basemap[:, 0], obs_pos_geo_basemap[:, 1])

    # # ###################################################################
    # ### discretize the fault segments
    # ### create slip basis vectors for each patch
    # ### build regularization matrix
    # ###################################################################
    Ds = len(fault_list[0]["slip_basis"]);   # number of basis vectors for slip patch, assumed same for all slip patches
    patches = [];  # a growing list of fault patches.
    total_fault_slip_basis = np.zeros((0, Ds, 3));   # a collection of basis vectors, one object for each fault patch
    patches_f = [];  # patches repeated for each basis (for example, one for dip slip and one for strike slip)
    slip_basis_f = np.zeros((0, 3))   # basis functions repeated for each slip patch
    fault_names_array = [];  # a list of fault names (integers) for each fault patch
    L_array = [];   # may hold several smoothing matrices, if using 2+ faults

    # # Fault processing
    for fault in fault_list:
        # Convert fault to cartesian coordinates
        fault["seg_pos_cart"] = plotting_library.geodetic_to_cartesian(fault["seg_pos_geo"], bm)

        # Discretize this fault segment
        seg = slippy.patch.Patch(fault["seg_pos_cart"], fault["length"], fault["width"], fault["strike"], fault["dip"]);
        single_fault_patches = np.array(seg.discretize(fault["Nlength"], fault["Nwidth"]))
        Ns = len(single_fault_patches);  # number of segments discretized from this fault

        # Create slip basis vectors, into an object of basis vectors for each fault segment
        single_fault_slip_basis = np.array([fault["slip_basis"] for _j in range(Ns)])
        total_fault_slip_basis = np.concatenate((total_fault_slip_basis, single_fault_slip_basis), axis=0)

        # Reshape: Packaging of slip_basis_f and patches_f
        single_fault_silp_basis_f = single_fault_slip_basis.reshape((Ns * Ds, 3))
        single_fault_patches_f = single_fault_patches[:, None].repeat(Ds, axis=1).reshape((Ns * Ds,))
        patches = np.concatenate((patches, single_fault_patches), axis=0)
        patches_f = np.concatenate((patches_f, single_fault_patches_f), axis=0)
        slip_basis_f = np.concatenate((slip_basis_f, single_fault_silp_basis_f), axis=0);
        names_for_patch = np.array([fault["name"]]).repeat(Ns);   # patch names: an integer repeated Ns times
        fault_names_array = np.concatenate((fault_names_array, names_for_patch), axis=0);

        ### build regularization matrix for this fault (for smoothing penalty?)
        L = np.zeros((0, Ns * Ds))
        indices = np.arange(Ns * Ds).reshape((Ns, Ds))
        for i in range(Ds):  # for each dimension of the basis
            connectivity = indices[:, i].reshape((fault["Nlength"], fault["Nwidth"]))
            Li = slippy.tikhonov.tikhonov_matrix(connectivity, 2, column_no=Ns * Ds)
            L = np.vstack((Li, L))

        L *= fault["penalty"]   # multiplying by smoothing strength for this fault
        L_array.append(L)   # collecting full smoothing matrix for each fault

    L = scipy.linalg.block_diag(*L_array)  # For 2+ faults: Make block diagonal matrix for tikhonov regularization

    # PARSE HOW MANY EPOCHS WE ARE USING
    # Tell us how many epochs, model parameters total, and output files we need.
    n_epochs = 0;
    total_spans, span_output_files = [], [];
    for epoch in config["epochs"].keys():
        n_epochs = n_epochs + 1;
        total_spans.append(config["epochs"][epoch]["name"]);
        span_output_files.append(config["output_dir"]+config["epochs"][epoch]["slip_output_file"]);
    n_model_params = np.shape(L)[0];    # model parameters that aren't leveling offset
    n_cols_bigG = n_model_params * n_epochs;  # the total number of fault-related model parameters across all time
    print("Finding fault model for: %d epochs " % n_epochs);
    print("Number of fault model parameters per epoch: %d" % n_model_params);
    print("Total number of model parameters: %d" % (n_model_params * n_epochs + number_of_datasets));  # seems singular?

    # INITIALIZING THE LARGE MATRICES
    d_total = np.zeros((0,));         # data vector
    sig_total = np.zeros((0,));       # uncertainties vector
    weight_total = np.zeros((0,));    # weights vector
    G_total = np.zeros((0, n_epochs * n_model_params));  # does not contain leveling offsets
    G_nosmooth_total = np.zeros((0, n_epochs * n_model_params));  # for computing predicted displacements
    nums_obs = [];    # list that holds number of data points in each dataset (3ngps, ninsar, etc.)
    pos_obs_list = [];     # list of [llh] for each observation, (ninsar + nlev + 3ngps)
    pos_basis_list = [];    # lists of look vectors or gps basis, (ninsar + nlev + 3ngps)

    # INITIAL DATA SCOPING: HOW MANY FILES WILL NEED TO BE READ?
    input_file_list, output_file_list = [], [];
    spans_list, strengths_list, signs_list = [], [], [];
    data_type_list, row_span_list = [], [];
    print("Available data indicated in json file: ")

    # Unpacking metadata for later use
    for data_file in config["data_files"].keys():
        input_file_list.append(config["data_files"][data_file]["data_file"]);  # 'infile' expected of all data files
        output_file_list.append(config["output_dir"]+config["data_files"][data_file]["outfile"]);  # predicted outfile
        data_type_list.append(config["data_files"][data_file]["type"]);       # 'type' expected of all data files
        spans_list.append(config["data_files"][data_file]["span"]);           # 'span' expected of all data files
        strengths_list.append(config["data_files"][data_file]["strength"]);   # 'strength' expected of all data files
        print(config["data_files"][data_file]["data_file"], " spans ", config["data_files"][data_file]["span"]);
        if config["data_files"][data_file]["type"] == "leveling":
            signs_list.append(config["data_files"][data_file]["sign"]);
        else:
            signs_list.append(0.0);   # june 2021: does this result in a singular matrix? It might...
            # we seem to have one parameter for each dataset, regardless of whether it needs it.

    # START THE MAJOR INVERSION LOOP
    for filenum in range(len(input_file_list)):

        # INPUT STAGE
        filename = input_file_list[filenum];
        if data_type_list[filenum] == 'gps':
            [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ngps] = input_gps_file(filename);
            obs_weighting_f = (1 / strengths_list[filenum]) * np.ones((Ngps * 3,));
        elif data_type_list[filenum] == 'insar' or data_type_list[filenum] == "leveling":
            [obs_disp_f, obs_sigma_f, obs_basis_f, obs_pos_geo_f, Ninsar] = input_insar_file(filename);
            obs_weighting_f = (1 / strengths_list[filenum]) * np.ones((Ninsar,));
        else:
            print("ERROR! Unrecognized data type %s " % data_type_list[filenum]);
            continue;

        # Metadata collection about this input source
        pos_obs_list.append(obs_pos_geo_f);  # lists of llh, for writing outputs later
        pos_basis_list.append(obs_basis_f);   # lists of look vectors or gps basis, for writing outputs later
        nums_obs.append(len(obs_disp_f));    # list append: number of obs in this dataset (ninsar, 3ngps, etc)
        this_data_spans = spans_list[filenum];
        print("This data spans:", this_data_spans);

        # # Geodetic to Cartesian coordinates for this dataset
        obs_pos_cart_f = plotting_library.geodetic_to_cartesian(obs_pos_geo_f, bm)

        # BUILD SYSTEM MATRIX
        # here, leveling=False because we'll add leveling manually later
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

        # BUILD THE G MATRIX FOR THIS SET OF OBSERVATIONS
        # (L DEPENDS ON FAULT GEOMETRY ONLY)
        G_ext, d_ext = G_with_smoothing(G, L, alpha, obs_disp_f);
        print("Shape of G, L:", np.shape(G), " ", np.shape(L))
        print("Shape of Gext (G,L,alpha):", np.shape(G_ext));

        # APPEND TO THE BIGGER MATRIX AND DATA VECTOR (WITH THE CORRECT SPANS)
        G_rowblock_obs = np.zeros((np.shape(G_ext)[0], n_cols_bigG));  # one block per data file
        G_nosmoothing_obs = np.zeros((np.shape(G)[0], n_cols_bigG));  # one block per data file, for forward modeling

        # Which spans does the data cover? This is like the "1" in SBAS
        # One new "rowblock" for each data file (blocks that cross all the time spans)
        count = 0;
        for epoch in total_spans:
            col_limits = [count * n_model_params, (count + 1) * n_model_params];
            if epoch in this_data_spans:
                G_rowblock_obs[:, col_limits[0]:col_limits[1]] = G_ext;   # fill rowblock with copies of G_ext block
                G_nosmoothing_obs[:, col_limits[0]:col_limits[1]] = G;  # fill nosmooth-rowblock with copies of G block
            count = count + 1;

        # APPEND NEW ROWBLOCK TO G_TOTAL AND G_NOSMOOTH_TOTAL
        d_total = np.concatenate((d_total, d_ext));  # building up the total data vector
        sig_total = np.concatenate((sig_total, obs_sigma_f));  # building up the total sigma vector
        weight_total = np.concatenate((weight_total, obs_weighting_f));  # building up the total weighting vector
        top_row_data = len(G_total);  # Where in the matrix rows is this data?
        bottom_row_data = len(G_total) + len(G);  # Where in the matrix rows is this data?
        row_span_list.append([top_row_data, bottom_row_data]);
        G_total = np.concatenate((G_total, G_rowblock_obs));
        print("  Adding %d lines " % len(G_rowblock_obs))
        G_nosmooth_total = np.concatenate((G_nosmooth_total, G_nosmoothing_obs));

    # End Read_Data and Build_G stage

    # ADDING THE COLUMNS FOR LEVELING OFFSETS TO G_TOTAL MATRIX (Corresponding to the data lines only)
    print("------\nBefore adding lines for leveling offsets, shape(G): ", np.shape(G_total));
    numrows = np.shape(G_total)[0];
    num_rows_nosmooth = np.shape(G_nosmooth_total)[0];
    count = 0;
    for filenum in range(len(input_file_list)):
        leveling_sign = signs_list[filenum];
        newcol = np.zeros((numrows, 1));
        for i in range(len(newcol)):
            if row_span_list[filenum][0] < i < row_span_list[filenum][1]:
                newcol[i] = leveling_sign;
        G_total = np.hstack((G_total, newcol));
        print("Adding column for %s " % input_file_list[filenum])

        # ADDING THE COLUMNS TO G_NOSMOOTH_TOTAL
        newcol_nosmooth = np.zeros((num_rows_nosmooth, 1));
        for i in range(len(newcol_nosmooth)):
            if count <= i < count + (row_span_list[filenum][1] - row_span_list[filenum][0]):
                newcol_nosmooth[i] = leveling_sign;
        count = count + (row_span_list[filenum][1] - row_span_list[filenum][0]);
        G_nosmooth_total = np.hstack((G_nosmooth_total, newcol_nosmooth));

    # Running a resolution test if desired.
    if int(config["resolution_test"]):
        resolution_tests.analyze_model_resolution_matrix(G_nosmooth_total, config["output_dir"]);

    # INVERT BIG-G
    #   ### estimate slip and compute predicted displacement
    #   #####################################################################
    slip_f = reg_nnls(G_total, d_total)  # the model
    print("G_total:", np.shape(G_total));
    print("slip_f:", np.shape(slip_f))
    print("G_nosmooth_total:", np.shape(G_nosmooth_total));
    print("shape of data (G*slip): ", np.shape(G_nosmooth_total.dot(slip_f)))
    print("shape of sigmas       : ", np.shape(sig_total))
    print("shape of weights      : ", np.shape(weight_total));
    pred_disp_f = G_nosmooth_total.dot(slip_f) * sig_total * weight_total;

    Ns_total = len(patches);
    total_cardinal_slip, leveling_offsets = parse_slip_outputs(slip_f, Ns_total, Ds, n_epochs, total_fault_slip_basis,
                                                               number_of_datasets);
    disp_models = parse_disp_outputs(pred_disp_f, nums_obs);

    # Defensive programming (Reporting errors)
    for i in range(len(leveling_offsets)):
        if signs_list[i] != 0:
            print("Leveling Offset for %s = %f m " % (input_file_list[i], leveling_offsets[i]));
            if abs(leveling_offsets[i]) < 0.0000001:
                print("WARNING: Leveling offset for %s close to zero. Consider a negative offset in G." %
                      (input_file_list[i]));

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
        slip_output_file = span_output_files[i];

        # ### write output
        # #####################################################################
        slippy.io.write_slip_data(patches_pos_geo,
                                  patches_strike, patches_dip,
                                  patches_length, patches_width,
                                  total_cardinal_slip[i], fault_names_array, slip_output_file)
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
            pred_disp_insar = disp_models[filenum];
            slippy.io.write_insar_data(pos_obs_list[filenum],
                                       pred_disp_insar, 0.0 * pred_disp_insar,
                                       pos_basis_list[filenum],
                                       output_file_list[filenum])
            print("Writing file %s " % output_file_list[filenum]);

        elif data_type_list[filenum] == 'leveling':
            pred_disp_leveling = disp_models[filenum];
            slippy.io.write_insar_data(pos_obs_list[filenum],
                                       pred_disp_leveling, 0.0 * pred_disp_leveling,
                                       pos_basis_list[filenum],
                                       output_file_list[filenum])
            print("Writing file %s " % output_file_list[filenum]);

    # MISC OUTPUTS: Graph of big G (with all smoothing parameters inside)
    graph_big_G(config, G_total);
    return;
