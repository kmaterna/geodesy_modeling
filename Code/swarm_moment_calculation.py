
# The purpose of these tools is to help read earthquake catalogs and do moment calculations

import numpy as np
import matplotlib.pyplot as plt  
import collections
import datetime as dt
import moment_calculations

Catalog = collections.namedtuple('Catalog',['lon','lat','depth','mag','fm']);

def read_Wei_2015_supplement(filename):
	print("Reading earthquake catalog from file %s " % filename)
	lon=[]; lat=[]; depth=[]; mag=[]; fm=[];
	ifile=open(filename);
	for line in ifile:
		lon.append(float(line.split()[2]));
		lat.append(float(line.split()[1]));
		depth.append(float(line.split()[3]));
		mag.append(float(line.split()[4]));
		fm.append(line.split()[5]);
	ifile.close();
	MyCatalog = Catalog(lon=lon, lat=lat, depth=depth, mag=mag, fm=fm);
	return MyCatalog;


def wei_catalog_total_moment(filename):
	# Shengji's 2015 supplementary material had a nice catalog of the major events in the 2012 Brawley Sequence. 
	# What was the total seismic moment released?  
	# Looks like it was about equivalent to a 5.57
	# This is significantly lower than what I invert geodetically for the swarm. 
	MyCat = read_Wei_2015_supplement(filename);
	M_total = 0;
	for i in range(len(MyCat.lon)):
		Mo = moment_calculations.moment_from_mw(MyCat.mag[i]);
		print("Adding Mw %f " % MyCat.mag[i]);
		M_total = M_total + Mo;
	print("Total moment: %f N-m" % M_total);
	Mw_total = moment_calculations.mw_from_moment(M_total);
	print("Total moment is %f N-m, corresponding to a Mw %f" % (M_total, Mw_total) );
	return; 
