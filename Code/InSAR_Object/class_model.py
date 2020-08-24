import collections

# FORMAT DEFINITION
InSAR_Object = collections.namedtuple('InSAR_Object', ['lon', 'lat', 'LOS', 'LOS_unc',
                                                       'lkv_E', 'lkv_N', 'lkv_U', 'starttime', 'endtime']);
# A generalized InSAR format
# displacements in mm
# lon, lat, LOS, and LOS_unc are vectors
# Look vector is which direction? Ground to satellite I believe.
# starttime and endtime are just for metadata
