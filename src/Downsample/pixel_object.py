# the downsampled pixel object

class Downsampled_pixel:
    """
    An object containing a quadtree-downsampled pixel, including
    downsampled pixel footprint, downsampled pixel look vector, and downsampled deformation values.
    """
    def __init__(self, mean, median, std, BL_corner, TR_corner, unitE, unitN, unitU):
        self.mean = mean;  # mean of LOS values within the pixel (meters)
        self.median = median;  # median of LOS values within the pixel (meters)
        self.std = std;  # standard deviation of LOS values within the pixel (meters)
        self.BL_corner = BL_corner;  # Coordinates of Bottom Left corner (longitude, latitude)
        self.TR_corner = TR_corner;  # Coordinates of Top Right corner (longitude, latitude)
        self.unitE = unitE;  # east component of unit vector from ground to satellite
        self.unitN = unitN;  # north component of unit vector from ground to satellite
        self.unitU = unitU;  # up component of unit vector from ground to satellite

    def set_mean(self, new_mean):
        self.mean = new_mean;

    def set_median(self, new_median):
        self.median = new_median;

    def set_std(self, new_std):
        self.std = new_std;
