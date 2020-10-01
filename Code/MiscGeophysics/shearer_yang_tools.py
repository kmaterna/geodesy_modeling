

import numpy as np 
import datetime as dt 
import csv
import qtm_tools

def input_shearer_cat(filename):
	print("Reading file %s " % filename);
	ifile=open(filename);
	dtarray = []; latitude = []; longitude = []; depth = []; magnitude = [];
	for line in ifile:
		temp = line.split();
		year = temp[0]
		month = temp[1]
		day = temp[2]
		hour = temp[3]
		minute = temp[4] 
		second = temp[5][0:2]
		if int(second)>59:
			print("We found a leap second at: ",year, month, day, hour, minute, second);
			print("Turning it back one second")
			second = '59'
		if hour == '-1':
			print("We found a problem at: ",year, month, day, hour, minute, second);
			print("Turning it forward one second")
			hour = '00'
			minute = '00'
			second = '00'
		eqdate = dt.datetime.strptime(year+" "+month+" "+day+" "+hour+" "+minute+" "+second,"%Y %m %d %H %M %S");
		dtarray.append(eqdate);
		latitude.append(float(temp[7]))
		longitude.append(float(temp[8]))
		depth.append(float(temp[9]))
		magnitude.append(float(temp[10]));
	ifile.close();
	
	return MyCat;


def read_usgs_website_csv(filename):
	# When you hit the 'DOWNLOAD' button on the USGS website
	dtarray = []; latitude = []; longitude = []; depth = []; magnitude = [];
	with open(filename) as csvfile:
		mycatreader = csv.reader(csvfile);
		for row in mycatreader:
			if row[0]=='time':
				continue;
			dtarray.append(dt.datetime.strptime(row[0][0:19],"%Y-%m-%dT%H:%M:%S"));
			latitude.append(float(row[1]));
			longitude.append(float(row[2]));
			depth.append(float(row[3]));
			magnitude.append(float(row[4]));

	MyCat = qtm_tools.Catalog(dtarray=dtarray, lon=longitude, lat=latitude, depth=depth, Mag=magnitude);

	return MyCat

