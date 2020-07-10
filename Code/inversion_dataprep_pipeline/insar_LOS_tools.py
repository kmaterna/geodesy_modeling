# July 2020
# Perform uniform downsampling on an InSAR dataset
# Impose Bounding Box
# Convert between InSAR data formats
# Remove a ramp empirically

import numpy as np 
import sys
import collections

InSAR_Object = collections.namedtuple('InSAR_Object',['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime']);



def uniform_downsampling(InSAR_obj, sampling_interval, averaging_window):
	# InSAR Object is similar to the Hines format

	return 0;



def InSAR_bounding_box(InSAR_obj, bbox=[-180, 180, -90, 90]):
	# Impose a bounding box on some InSAR data
	lon=[]; lat=[]; LOS=[]; LOS_unc=[]; unit_E=[]; unit_N=[]; unit_U=[];
	for i in range(len(InSAR_obj.lon)):
		if InSAR_obj.lon[i]>=bbox[0] and InSAR_obj.lon[i]<=bbox[1]:
			if InSAR_obj.lat[i]>=bbox[2] and InSAR_obj.lat[i]<=bbox[3]:
				if np.isnan(InSAR_obj.LOS[i]):
					continue;
				else:
					lon.append(InSAR_obj.lon[i]);
					lat.append(InSAR_obj.lat[i]);
					LOS.append(InSAR_obj.LOS[i]);
					LOS_unc.append(InSAR_obj.LOS_unc[i]);
					unit_E.append(InSAR_obj.lkv_E[i]);
					unit_N.append(InSAR_obj.lkv_N[i]);
					unit_U.append(InSAR_obj.lkv_U[i]);
	newInSAR_obj = InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N, lkv_U=unit_U, 
		starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
	return newInSAR_obj;



def TRE_to_InSAR_Obj(TRE_obj):
	# Divide by the number of years of observation
	# Return the Vert and East as separate objects
	tdelta = TRE_obj.endtime-TRE_obj.starttime;
	tre_interval_years = tdelta.days / 365.24;  # the number of years spanned by the TRE velocity. 
	Vert_LOS = TRE_obj.vvel*tre_interval_years;
	East_LOS = TRE_obj.evel*tre_interval_years;
	zeros = np.zeros(np.shape(TRE_obj.vvel));
	ones = np.zeros(np.shape(TRE_obj.vvel));
	Vert_obj = InSAR_Object(TRE_obj=lon, TRE_obj=lat, LOS=Vert_LOS, LOS_unc=TRE_obj.vvel_std, lkv_E=zeros, lkv_N=zeros, lkv_U=ones, 
		starttime=TRE_obj.starttime, endtime=TRE_obj.endtime);
	East_obj = InSAR_Object(TRE_obj=lon, TRE_obj=lat, LOS=East_LOS, LOS_unc=TRE_obj.evel_std, lkv_E=ones, lkv_N=zeros, lkv_U=zeros, 
		starttime=TRE_obj.starttime, endtime=TRE_obj.endtime);
	# Return two objects here
	return Vert_obj, East_obj;

