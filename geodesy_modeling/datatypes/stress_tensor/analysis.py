
import numpy as np
import pandas as pd
from Tectonic_Utils.geodesy import fault_vector_functions as fvf
from elastic_stresses_py.PyCoulomb import conversion_math

"""
The object of the stress tensor that we operate upon is just a 3x3 numpy matrix.
No need to create a special class for this object. 
"""


def resolve_stress_geometries(tau, friction, B=0, increment=10):
    """
    Convert the six-component stress tensor into its shear, normal, and Coulomb traction components
    on every possible geometry. Return a Pandas array.

    :param tau: 3x3 matrix of stress tensor, in pa
    :param friction: float, coefficient of friction
    :param B: Skempton's coefficient
    :param increment: integer, increment that we step through strike, dip, rake
    :return: Pandas dataframe with six columns [strike, dip, rake, shear, normal, coulomb] in kpa
    """
    # Step 1: Resolve the stress tensor on geometries
    # initialize list
    stresstemp = []
    # loop over strike, dip and rake
    for strike in np.arange(0, 360, increment):
        for dip in np.arange(0, 100, increment):
            for rake in np.arange(-180, 180, increment):
                # First compute the geometric vectors associated with the receiver
                strike_unit_vector = fvf.get_strike_vector(strike)  # a 3d vector in the horizontal plane.
                dip_unit_vector = fvf.get_dip_vector(strike, dip)  # a 3d vector.
                plane_normal = fvf.get_plane_normal(strike, dip)  # a 3d vector.

                effective_normal_stress, shear_stress, coulomb_stress = \
                    conversion_math.get_coulomb_stresses_internal(tau, strike_unit_vector, rake, dip_unit_vector,
                                                                  plane_normal, friction, B)
                # Nex
                stresstemp.append([strike, dip, rake, shear_stress, effective_normal_stress, coulomb_stress])

    # turn results to dataframe
    stresses = pd.DataFrame(stresstemp, columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb'])
    return stresses
