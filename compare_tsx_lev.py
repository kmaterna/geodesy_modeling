# Compare TSX from 2012-2013 with Leveling
# Holy Cow I reproduced Mariana's results. 

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import datetime as dt 
import sys
import multiSAR_compare_tools
import multiSAR_input_functions
import haversine


def plot_TSX(tsxData):
	plt.figure(figsize=(10,10));
	plt.scatter(tsxData.lon, tsxData.lat, c=tsxData.vvel, s=10);
	plt.colorbar();
	plt.savefig("Comparisons/TSX/TSX_vvel.png");
	return;


def one_to_one_comparison(myLev, tsxData, vector_index):
	lev1=3;
	lev2=4;  # the TSX interval. 
	filename = "Comparisons/TSX/one_to_one_"+str(lev1)+str(lev2);

	oto_lev=[];
	oto_tsx=[];
	lon_plotting=[];
	lat_plotting=[];


	for i in range(len(myLev.lon)):
		if vector_index[i]==np.nan:
			continue;
		else:
			leveling_offset=1000*(myLev.leveling[i][lev2]-myLev.leveling[i][lev1]) # negative sign convention
			tsx_offset = tsxData.vvel[vector_index[i]];
			if ~np.isnan(leveling_offset) and ~np.isnan(tsx_offset):
				oto_lev.append(leveling_offset);
				oto_tsx.append(tsx_offset);  # This would typically be converted into a displacement by multiplying by years. 
				lon_plotting.append(myLev.lon[i]);
				lat_plotting.append(myLev.lat[i]);

	vmin=-50; vmax=50;
	fig,axarr = plt.subplots(2,2, figsize=(14,10));

	# Individual plot of leveling or TSX
	color_boundary_object=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object,cmap='RdYlBu_r');
	for i in range(len(oto_lev)):
		dot_color=custom_cmap.to_rgba(oto_lev[i]);
		dot_color_tsx=custom_cmap.to_rgba(oto_tsx[i]);
		axarr[0][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color, fillstyle="full");
		axarr[1][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color_tsx, fillstyle="full");

	axarr[0][0].set_title("Leveling: "+dt.datetime.strftime(myLev.dtarray[lev1],"%m-%Y")+" to "+dt.datetime.strftime(myLev.dtarray[lev2],"%m-%Y"),fontsize=15);
	axarr[0][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12,color='black');
	axarr[1][0].set_title("TSX: "+dt.datetime.strftime(tsxData.starttime,"%m-%Y")+" to "+dt.datetime.strftime(tsxData.endtime,"%m-%Y"),fontsize=15);
	axarr[1][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12,color='black');

	# The one-to-one plot
	axarr[0][1].plot([-80,80],[-80,80],linestyle='--',color='gray')
	axarr[0][1].plot(oto_lev, oto_tsx, markersize=10, marker='v',linewidth=0);
	axarr[0][1].set_xlabel('Leveling offset (mm)',fontsize=20);
	axarr[0][1].set_ylabel('TSX offset (mm)',fontsize=20);
	axarr[0][1].tick_params(axis='both',which='major',labelsize=16)
	axarr[0][1].set_xlim([-80,80])
	axarr[0][1].set_ylim([-80,80])
	axarr[0][1].grid(True)

	axarr[1][1].scatter(tsxData.lon, tsxData.lat, c=tsxData.vvel, s=8,
		marker='o',cmap='RdYlBu_r',vmin=vmin, vmax=vmax);
	axarr[1][1].plot(myLev.lon, myLev.lat, '*',color='black');
	axarr[1][1].plot(myLev.lon[0], myLev.lat[0], '*', color='red');
	axarr[1][1].plot(-115.510,33.081,'v',markersize=10,color='black');

	cbarax = fig.add_axes([0.75,0.35,0.2,0.3],visible=False);
	color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
	custom_cmap.set_array(np.arange(vmin, vmax));
	cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
	cb.set_label('Displacement (mm)', fontsize=18);
	cb.ax.tick_params(labelsize=12);

	plt.savefig(filename+".png");

	return;



if __name__=="__main__":
	file_dict = multiSAR_input_functions.get_file_dictionary();
	tsxData = multiSAR_input_functions.inputs_tsx(file_dict["tsx"]);
	myLev = multiSAR_input_functions.inputs_leveling(file_dict["leveling"].split()[0], file_dict["leveling"].split()[1]);
	myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);
	vector_index = multiSAR_compare_tools.find_leveling_in_vector(myLev, tsxData);

	one_to_one_comparison(myLev, tsxData, vector_index);
	plot_TSX(tsxData);



