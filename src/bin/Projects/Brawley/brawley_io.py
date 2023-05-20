

def get_file_dictionary(config_filename):
    """GET FILE NAMES from Brawley file dictionary."""
    this_dict = {};
    print("Reading file %s " % config_filename);
    ifile = open(config_filename);
    for line in ifile:
        data_type = line.split(':')[0];
        total_data_files = line.split()[1];  # assuming one file per list entry
        this_dict[data_type] = total_data_files;
    ifile.close();
    return this_dict;
