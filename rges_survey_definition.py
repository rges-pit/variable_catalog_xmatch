# Utility functions providing information describing the
# Roman Galactic Exoplanet Survey of the Galactic Bulge

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

class RGESSurvey():
    """
    Class describing the paraameters of the Roman Galactic Exoplanet Survey
    """
    def __init__(self):
        self.field_radius = 4.0 * u.deg
        self.field_centers = [
        (267.835895375*u.deg, -30.0608178195*u.deg)
        ]
        self.fields = [SkyCoord(x[0], x[1], frame='icrs', unit=(u.deg, u.deg)) for x in self.field_centers]

    def chk_in_field(self, ifield, s):
        """
        Method to verify if a given RA, Dec lies in an individual pointing
        footprint

        :param ifield: Index of field in fields list
        :param s: SkyCoord of the RA, Dec to check
        :return: Boolean
        """

        f = self.fields[ifield]
        if f.separation(s) <= self.field_radius:
            return True
        else:
            return False

    def chk_in_survey(self, ra, dec):
        """
        Method to verify if a given RA, Dec lies within the
        RGES survey footprint

        :param ra float, deg Right Ascension
        :param dec float, deg Declination
        :returns Boolean
        """

        s = SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))

        field_checks = [self.chk_in_field(i,s) for i in range(0,len(self.fields),1)]

        return any(field_checks)

    def find_stars_in_survey(self, coord_list):
        """
        Method to identify stars from the coordinate list given that lie within
        the RGES survey footprint

        :param coord_list: Multi-object SkyCoord
        :return: array of indices of objects within the footprint
        """

        survey_stars = np.array([], dtype=int)
        for f in self.fields:
            sep = f.separation(coord_list)
            idx = np.where(sep <= self.field_radius)[0]
            if len(idx) > 0:
                survey_stars = np.concatenate((survey_stars, idx))

        return survey_stars.tolist()