import collections

# FORMAT DEFINITION
InSAR_2D_Object = collections.namedtuple('InSAR_2D_Object', ['lon', 'lat', 'LOS', 'LOS_unc', 'lkv_E', 'lkv_N',
                                                             'lkv_U', 'starttime', 'endtime']);

"""
A generalized grid InSAR format
Displacements in mm (if LOS is a displacement measurement instead of phase or other)
Look vector is Ground to Satellite.
starttime and endtime are just for metadata
lon, lat are 1D arrays
LOS, LOS_unc are a 2D arrays
lkv_E, lkv_N, lkv_U are 2D arrays. 
"""
