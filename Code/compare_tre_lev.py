# Compare TRE Data with Leveling: 
# TSX from 2012-2013
# S1 from 2014-2018
# S1 from 2018-2019
# Holy Cow I reproduced Mariana's results. 

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import datetime as dt 
import sys
import multiSAR_utilities
import multiSAR_input_functions
import InSAR_Object.inputs


def plot_TRE(tsxData,output_dir):
	plt.figure(figsize=(10,10));
	plt.scatter(tsxData.lon, tsxData.lat, c=tsxData.vvel, s=10);
	plt.colorbar();
	plt.savefig(output_dir+"TRE_vvel.png");
	return;


def one_to_one_comparison(myLev, treData, vector_index, lev1, lev2, sat, output_dir):
	filename = output_dir+"one_to_one_"+str(lev1)+str(lev2);

	oto_lev=[];
	oto_tsx=[];
	lon_plotting=[];
	lat_plotting=[];
	tdelta = treData.endtime-treData.starttime;
	tre_interval = tdelta.days / 365.24;  # the number of years spanned by the TRE velocity. 


	for i in range(len(myLev.lon)):
		if vector_index[i]==np.nan:
			continue;
		else:
			leveling_offset=1000*(myLev.leveling[i][lev2]-myLev.leveling[i][lev1]) # negative sign convention
			tre_offset = treData.vvel[vector_index[i]] * tre_interval;
			if ~np.isnan(leveling_offset) and ~np.isnan(tre_offset):
				oto_lev.append(leveling_offset);
				oto_tsx.append(tre_offset);  # This would typically be converted into a displacement by multiplying by years. 
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
	axarr[1][0].set_title(sat+": "+dt.datetime.strftime(treData.starttime,"%m-%Y")+" to "+dt.datetime.strftime(treData.endtime,"%m-%Y"),fontsize=15);
	axarr[1][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12,color='black');

	# The one-to-one plot
	axarr[0][1].plot([-80,80],[-80,80],linestyle='--',color='gray')
	axarr[0][1].plot(oto_lev, oto_tsx, markersize=10, marker='v',linewidth=0);
	axarr[0][1].set_xlabel('Leveling offset (mm)',fontsize=20);
	axarr[0][1].set_ylabel(sat+' offset (mm)',fontsize=20);
	axarr[0][1].tick_params(axis='both',which='major',labelsize=16)
	axarr[0][1].set_xlim([-80,80])
	axarr[0][1].set_ylim([-80,80])
	axarr[0][1].grid(True)

	tre_disps = [i*tre_interval for i in treData.vvel];
	axarr[1][1].scatter(treData.lon, treData.lat, c=tre_disps, s=8,
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
	# Opening stuff
	config_filename = "config_file.txt"
	file_dict = multiSAR_input_functions.get_file_dictionary(config_filename);
	myLev = multiSAR_input_functions.inputs_leveling(file_dict["leveling"].split()[0], file_dict["leveling"].split()[1]);
	myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);

	#TSX experiment
	output_dir = "TSX/"
	tsxData = InSAR_Object.inputs.inputs_TRE(file_dict["tsx"]);
	vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, tsxData);
	one_to_one_comparison(myLev, tsxData, vector_index, 3, 4, "TSX", output_dir);  # the 2012-2013 interval
	plot_TRE(tsxData,output_dir);

	# S1 experiment
	output_dir = "SNT1/"
	s1Data = InSAR_Object.inputs.inputs_TRE(file_dict["snt1"]);
	vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, s1Data);
	one_to_one_comparison(myLev, s1Data, vector_index, 5, 8, "S1", output_dir); # the 2014-2018 interval
	plot_TRE(s1Data,output_dir);

	# S1 experiment
	output_dir = "SNT2/"
	s1Data = InSAR_Object.inputs.inputs_TRE(file_dict["snt2"]);
	vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, s1Data);
	one_to_one_comparison(myLev, s1Data, vector_index, 8, 9, "S1", output_dir); # the 2014-2018 interval
	plot_TRE(s1Data,output_dir);

