# Test suite for the retrieval and handling of OGLE lightcurve data

import pytest
from astropy.table import Table, Column
import utils
import numpy as np

@pytest.mark.parametrize(
    "test_star, expected",
    [
        (
        Table([
            Column(name='Name', data=['OGLE-BLG-LPV-023363']),
            Column(name='Type', data=['lpv']),
            Column(name='RA', data=[265.9569166666666]),
            Column(name='Dec', data=[-33.408722222222224])
        ]),
        [
            Table([
            Column(name='HJD', data=np.linspace(5656.58653, 7656.58653, 1000)),
            Column(name='mag', data=np.array([15.5]*1000)),
            Column(name='mag_error', data=np.array([0.005]*1000)),
        ]),
            None
        ]
        )
    ]
)
def test_fetch_ogle_photometry(test_star, expected):
    ilc, vlc = utils.fetch_ogle_photometry(test_star)

    # The test star given has an I-band lightcurve but no V-band data,
    # so ilc returns a table but vlc returns None
    assert(type(ilc) == type(expected[0]))
    assert(type(vlc) == type(expected[1]))
    assert(ilc.colnames == expected[0].colnames)