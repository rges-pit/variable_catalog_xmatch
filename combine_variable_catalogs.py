# Program to produce a combined catalog of known variables within the
# Roman Galactic Exoplanet Survey footprint, crossmatching between catalogs to
# identify stars where data is available from multiple sources

from os import path
import json
import argparse
import utils
import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def crossmatch_catalogs(args):
    """
    Function to combine variable star catalogs from OGLE and VVV into a single
    catalog, crossmatching to identify any common entries

    :param args:
    :return: None, outputs combined JSON catalog
    """

    # Load input catalogs, ensuring the coordinates are in a consistent format
    ogle_catalog = utils.load_json_catalog_as_table(args.ogle_cat_file, decimal_degrees=True)
    vvv_catalog = utils.load_json_catalog_as_table(args.vvv_cat_file)
    gaia_catalog = utils.load_json_catalog_as_table(args.gaia_cat_file, decimal_degrees=True)

    # The OGLE catalog is likely to be the most comprehensive, so we start with this,
    # and crossmatch all other catalogs against it.  To combine the information from
    # the other catalogs, we add key-value pairs to include the identifiers and types from
    # the other surveys, if a star is known to more than one.  We first extend the
    # combined_catalog to add columns to ensure this information is held for all
    # objects
    combined_catalog = copy.deepcopy(ogle_catalog)
    combined_catalog.add_column([None]*len(combined_catalog), name='VVV_ID')
    combined_catalog.add_column([None]*len(combined_catalog), name='VVV_type')
    combined_catalog.add_column([None]*len(combined_catalog), name='Gaia_ID')
    combined_catalog.add_column([None]*len(combined_catalog), name='Gaia_type')
    combined_catalog.add_column([None]*len(combined_catalog), name='UKIRT_source_table')
    combined_catalog.add_column([None]*len(combined_catalog), name='UKIRT_lc_files')
    print('Initial catalog has ' + str(len(combined_catalog)) + ' entries from OGLE')

    # Cross-match VVV catalog, combining entries if the stars are known to OGLE,
    # or adding them to the combined catalog if not
    merge_catalog(combined_catalog, vvv_catalog, 'VVV', match_radius_arcsec=2.0)
    merge_catalog(combined_catalog, gaia_catalog, 'Gaia', match_radius_arcsec=2.0)

    # Output combined catalog
    utils.output_json_catalog_from_table(combined_catalog, args.output_file)

def merge_catalog(combined_catalog, catalog, catalog_id, match_radius_arcsec=2.0):
    """
    Function to crossmatch entries from the new survey's catalog with those in the
    existing combined_catalog, which is taken as the main reference.

    :param combined_catalog: Table object
    :param catalog: Table object
    :return: combined_catalog: Table object
    """
    match_radius = (match_radius_arcsec / 3600.0)*u.deg

    # Produce an array of SkyCoords for the whole reference catalog
    ref_catalog = SkyCoord(combined_catalog['RA'], combined_catalog['Dec'],
                           frame='icrs', unit=(u.deg, u.deg))

    # For each entry in the second catalog, check for matching entries in the
    # reference catalog.  If no matches are found within the allowed match radius,
    # add a new entry to the combined catalog
    for star in catalog:
        s = SkyCoord(star['RA'], star['Dec'], frame='icrs', unit=(u.deg, u.deg))
        sep = s.separation(ref_catalog)
        idx = np.where(sep <= match_radius)[0]
        if len(idx) > 0:
            if catalog_id == 'VVV':
                combined_catalog[idx[0]]['VVV_ID'] = star['Name']
                combined_catalog[idx[0]]['VVV_type'] = star['Type']
            elif catalog_id == 'Gaia':
                combined_catalog[idx[0]]['Gaia_ID'] = star['Name']
                combined_catalog[idx[0]]['Gaia_type'] = star['Type']
            else:
                raise IOError('Unrecognised catalog ID: ' + catalog_id)
        else:
            if catalog_id == 'VVV':
                # When adding a VVV star to the catalog, the main 'Name' and 'Type'
                # columns are set to their VVV values, and these are echoed in the
                # VVV columns to indicate where this information comes from
                combined_catalog.add_row([
                    star['Name'],
                    star['Type'],
                    star['RA'],
                    star['Dec'],
                    star['Name'],
                    star['Type'],
                    None,           # Null entries for Gaia identifiers
                    None,
                    None,           # Null entries for UKIRT cross-matching
                    None
                ])
            elif catalog_id == 'Gaia':
                # Similarly, when adding Gaia entries, we use the Name and Type columns
                # from the crossmatched catalog
                combined_catalog.add_row([
                    star['Name'],
                    star['Type'],
                    star['RA'],
                    star['Dec'],
                    None,           # Null entries for VVV identifiers
                    None,
                    star['Name'],
                    star['Type'],
                    None,  # Null entries for UKIRT cross-matching
                    None
                ])
            else:
                raise IOError('Unrecognised catalog ID: ' + catalog_id)

    print('Length of combined catalog after merge: ' + str(len(combined_catalog)))

    return combined_catalog

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('ogle_cat_file', help='Path to OGLE catalog file')
    parser.add_argument('vvv_cat_file', help='Path to VVV catalog file')
    parser.add_argument('gaia_cat_file', help='Path to Gaia catalog file')
    parser.add_argument('output_file', help='Path to output combined catalog file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    crossmatch_catalogs(args)