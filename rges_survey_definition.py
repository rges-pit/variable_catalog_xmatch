# Utility functions providing information describing the
# Roman Galactic Exoplanet Survey of the Galactic Bulge

from astropy import units as u

def get_rges_field_centers():
    """
    Function provides the set of nominal field centers of the Roman Galactic
    Exoplanet Survey of the Galactic Bulge

    Although the Roman WFI footprint is irregular for the time being, this
    function approximates the field of view as a circle of radius ~10arcmin.

    :return: rges_fields  list of tuples (RA center, Dec center, radius)
    """

    # Default field radius, degrees
    field_radius = 0.168 * u.deg

    rges_fields = [
        (267.835895375*u.deg, -30.0608178195*u.deg, field_radius)
        ]

    return rges_fields