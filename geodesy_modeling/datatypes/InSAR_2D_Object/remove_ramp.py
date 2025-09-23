from ..InSAR_1D_Object import remove_ramp
from .class_model import Insar2dObject


def fit_ramp(InSAR_Obj: Insar2dObject):
    """
    Fit the best-fitting plane parameters to an object of InSAR data (2d).

    :param InSAR_Obj: 2D InSAR object
    :return: a, b, c
    """
    data1d = InSAR_Obj.convert_to_insar1D()
    a, b, c = remove_ramp.fit_ramp(data1d)
    return a, b, c


def remove_best_fit_ramp(InSAR_Obj: Insar2dObject):
    """
    Find the best-fitting ramp from a set of InSAR observations and remove it from the data.
    Plane equation: ax + by + c = z
    Solving Ax = B

    :param InSAR_Obj: 2D insar object
    :returns: 2D insar object
    """
    a, b, c = fit_ramp(InSAR_Obj)  # Solve for the best-fitting ramp of equation ax + by + c = z
    new_insar = InSAR_Obj.subtract_ramp(a, b, c)  # Removing the planar model
    return new_insar
