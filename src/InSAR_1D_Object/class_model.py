import collections

# FORMAT DEFINITION
InSAR_1D_Object = collections.namedtuple('InSAR_1D_Object', ['lon', 'lat', 'LOS', 'LOS_unc',
                                                             'lkv_E', 'lkv_N', 'lkv_U', 'starttime', 'endtime']);
"""
A generalized 1D InSAR format
displacements in mm
lon, lat, LOS, and LOS_unc are vectors
Look vector is Ground to Satellite
starttime and endtime are just for metadata
lkv_E, lkv_N, lkv_U are vectors
"""
