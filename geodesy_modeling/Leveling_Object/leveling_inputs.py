"""
Several functions (research code) for io of leveling data, including CEC Salton Trough Leveling Data
Some input functions also exist in specific project directories instead of being consolidated here
"""

import pandas
import os
import datetime as dt
import numpy as np
from .LevStation import LevStation


def inputs_brawley_leveling(data_filename, errors_filename):
    """
    Read leveling from CEC Salton Trough North Brawley leveling data
    Yes this is research code, but I'll keep this inside in case the leveing data type changes.
    """
    print("Reading in %s" % data_filename)
    df = pandas.read_excel(data_filename, engine='openpyxl')
    column_names = df.columns[1:-1].tolist()

    # Fix typos in metadata and data
    [corrtype, _sheet, rownum, colnum, _old_val, new_values] = read_cec_leveling_errors(errors_filename, data_filename)
    df = implement_changes_dataframe(df, corrtype, rownum, colnum, new_values)
    column_names = implement_changes_colnames(column_names, corrtype, rownum, colnum, new_values)

    # Reading name information, Fix Data and Metadata typos
    dtarray = get_datetimes(column_names)
    names = df['BENCHMARK'].values.tolist()

    # Reading lat/lon information
    lonlat_sheet = pandas.read_excel(data_filename, engine='openpyxl', sheet_name=2)
    ll_names = lonlat_sheet['Benchmark'].values.tolist()
    longitudes = [float(x) for x in lonlat_sheet['Longitude'].values.tolist()]
    latitudes = [float(x) for x in lonlat_sheet['Latitude'].values.tolist()]
    names, lons, lats = match_lon_lat(names, latitudes, longitudes, ll_names)

    LevStationList = []
    for x in range(0, 83):
        single_leveling_array = df.loc[x][1:-1].tolist()
        single_leveling_array = clean_single_ts(single_leveling_array)
        new_object = LevStation(name=names[x], lon=lons[x], lat=lats[x], dtarray=dtarray,
                                leveling=single_leveling_array, reflon=lons[0], reflat=lats[0])
        LevStationList.append(new_object)
    return LevStationList


def read_cec_leveling_errors(error_filename, data_filename):
    """Read file that documents typos and errors"""
    data_filename = os.path.split(data_filename)[1]
    print("Reading documented errors in %s " % error_filename)
    corrtype, sheetnum, rownum, colnum, old_values, new_values = [], [], [], [], [], []
    ifile = open(error_filename, 'r')
    for line in ifile:
        temp = line.split("::")
        if temp[0] == "#":
            continue
        if temp[0] == data_filename:  # filter to errors that apply to a particular data file
            if temp[4] == "nan":
                continue  # ignore a case for notes (not a real correction)
            corrtype.append(temp[1])  # data or metadata (basically: string or float?)
            sheetnum.append(int(temp[2]))
            rownum.append(int(temp[3]))
            colnum.append(int(temp[4]))
            old_values.append(temp[5])
            new_values.append(temp[6])
    ifile.close()
    return [corrtype, sheetnum, rownum, colnum, old_values, new_values]


def implement_changes_dataframe(df, corrtype, rownum, colnum, new_values):
    """Implement Data changes to array of data"""
    print("Implementing changes to data:")
    for i in range(len(rownum)):
        if corrtype[i] == 'Metadata':
            continue
        col_names = df.columns[0:-1].tolist()
        print("Finding error in column %s" % col_names[colnum[i]-1])  # get the right column
        old_value = df.iloc[rownum[i]-2, colnum[i]-1]
        df.iloc[rownum[i]-2, colnum[i]-1] = new_values[i]  # assignment
        print("   Carefully replacing data %s with %s" % (old_value, new_values[i]))
    return df


def implement_changes_colnames(column_names, corrtype, rownum, colnum, new_values):
    """Implement metadata changes to column names"""
    print("Implementing changes to metadata:")
    for i in range(len(rownum)):
        if corrtype[i] == 'Metadata' and rownum[i] == 0:  # if using a column name
            print("   Swapping \"%s\" for \"%s\"" % (column_names[colnum[i]-1], new_values[i]))
            column_names[colnum[i]-1] = new_values[i]   # one-indexed column numbers I guess
    return column_names


def get_datetimes(timestrings):
    dtarray = []
    for i in range(len(timestrings)):
        # Normal dates
        if " 88 " in timestrings[i]:
            temp = timestrings[i].split(" 88 ")
            temp2 = temp[1].split()
            mmm = temp2[0]
            year = temp2[1]
            dtarray.append(dt.datetime.strptime(year + "-" + mmm + "-01", "%Y-%b-%d"))  # issue here, but not too bad.
        else:  # For the special cases
            if "NOLTE 2008" in timestrings[i]:
                dtarray.append(dt.datetime.strptime("2008-Nov-01", "%Y-%b-%d"))
    return dtarray


def match_lon_lat(names, lats, lons, ll_names):
    """Pair up the latlon info with the timeseries info"""
    matched_lons, matched_lats = [], []
    for i in range(len(names)):
        find_name = names[i]
        if names[i] == "Y-1225 Datum":
            find_name = "Y 1225"
        idx = ll_names.index(find_name)
        matched_lats.append(lats[idx])
        matched_lons.append(lons[idx])
    return [names, matched_lons, matched_lats]


def clean_single_ts(array):
    newarray = []
    for i in range(len(array)):
        if str(array[i]) == "-" or str(array[i]) == "DESTROYED" or str(array[i]) == "DAMAGED" or str(
                array[i]) == "NOT" or str(array[i]) == "FOUND":
            newarray.append(np.nan)
        else:
            newarray.append(float(array[i]))
    return newarray


# LEVELING COMPUTE FUNCITON (REFERENCE TO DATUM)
def compute_rel_to_datum_nov_2009(data):
    """Skips the 2008 measurement. Returns an object that is 83x10"""
    RefLevStations = []
    for station in data:

        # Automatically find the first day that matters.  Either after 2008 or has data.
        idx = 0
        for j in range(len(station.dtarray)):
            if ~np.isnan(station.leveling[j]) and station.dtarray[j] > dt.datetime.strptime("2009-01-01", "%Y-%m-%d"):
                idx = j  # this is the first date after 2009 that has data
                break

        # Accounting for a change in Datum height in 2014
        idx_early = 6  # the placement of 2014 before adjustment on the spreadsheet
        idx_late = 7  # the placement of 2014 after adjustment on the spreadsheet
        step = station.leveling[idx_early] - station.leveling[idx_late]

        referenced_dates, referenced_data = [], []
        for j in range(1, len(station.dtarray)):  # skipping 2008 anyway.
            if j == 6:
                continue  # passing over the 2014 measurement before re-referencing.
            if station.dtarray[j] > dt.datetime.strptime("2014-01-01", "%Y-%m-%d"):
                referenced_dates.append(station.dtarray[j])
                referenced_data.append(station.leveling[j] - station.leveling[idx] + step)
            else:
                referenced_dates.append(station.dtarray[j])
                referenced_data.append(station.leveling[j] - station.leveling[idx])

        referenced_object = LevStation(name=station.name, lon=station.lon, lat=station.lat, dtarray=referenced_dates,
                                       leveling=referenced_data, reflon=station.reflon, reflat=station.reflat)
        RefLevStations.append(referenced_object)
    return RefLevStations


# HEBER DATA SPREADSHEET
def inputs_leveling_heber(infile):
    """
    CEC HEBER LEVELING SPREADSHEET INTO LIST OF LEVELING OBJECTS.
    """
    station_list = []
    print("Reading in %s" % infile)

    # Get locations of benchmarks and reference benchmark, with the reference in the first line.
    df = pandas.read_excel(infile, sheet_name=1)
    locnames, all_lons, all_lats = [], [], []
    names = df['BENCHMARKS HEBER'].values.tolist()
    lats_temp = df['Lat'].values.tolist()
    lons_temp = df['Lon'].values.tolist()
    for i in range(len(lats_temp)):
        if lats_temp[i] != "":
            latstring = str(lats_temp[i]).replace('..', '.')
            lonstring = str(lons_temp[i]).replace('..', '.')
            locnames.append(names[i])
            all_lats.append(float(latstring))
            all_lons.append(float(lonstring))
    reflat = all_lats[0]
    reflon = all_lons[0]

    # # Get leveling data from the spreadsheet
    df = pandas.read_excel(infile, sheet_name=0)
    dtstrings = df.iloc[31:56, 0].tolist()
    dtarray = [dt.datetime.strptime(x, "%b %Y") for x in dtstrings]

    # Extract each station's leveling data
    for colnum in range(5, 163):  # for each station's leveling data
        station_name = df.iloc[30][colnum]  # the string station name
        levarray = []
        for i in range(31, 56):
            if df.iloc[i][colnum] in ["DESTROYED", "", "NOT FOUND", "?", "UNACCESSABLE ", "LOST"]:
                levarray.append(np.nan)
            else:
                levarray.append(float(df.iloc[i][colnum]))
        station_lon_idx = locnames.index(station_name)
        new_station = LevStation(name=station_name, lat=all_lats[station_lon_idx], lon=all_lons[station_lon_idx],
                                 dtarray=dtarray, leveling=levarray, reflon=reflon, reflat=reflat)
        station_list.append(new_station)
    print("Returning %d leveling stations " % len(station_list))
    return station_list
