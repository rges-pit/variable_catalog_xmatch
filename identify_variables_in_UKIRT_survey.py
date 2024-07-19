# Program to identify UKIRT lightcurves, if available, for variables identified by
# other Bulge surveys including OGLE, VVV
# While it is possible to search for individual stars in the UKIRT data through the
# NExSci archive website, bulk queries for large numbers of stars require this
# automated approach

import argparse
import utils
import numpy as np
from os import path
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob

def find_ukirt_data_for_variables(args):
    """
    Function to search the source tables of the UKIRT microlensing survey for entries
    corresponding to stars from an input catalog.

    Originally this function searched for each known variable over all source tables,
    but since the source tables are so large they cannot be easily stored in memory,
    this approach demands we repeatedly reload the source tables, which is very slow.
    So the function was refactored to scan each table once for all known variables.

    :param args: ArgumentParser object
    :return: None, outputs updated JSON catalog with lightcurve files where available
    """

    ## NOTE TO SELF
    ## Also need to write a directory walker to catalog the paths to each lightcurve file,
    ## since it doesn't relate to the source table

    # Load the input catalog of stars to search for.  Add column(s) to store
    # UKIRT cross-matching information
    var_catalog = utils.load_json_catalog(args.catalog)
    if 'UKIRT_source_table' not in var_catalog.colnames:
        var_catalog.add_column([None] * len(var_catalog), name='UKIRT_source_table')
        var_catalog.add_column([None] * len(var_catalog), name='UKIRT_lc_files')

    # Load UKIRT look-up table of field pointings
    lut = utils.load_ukirt_lut(args.lut)

    # Make list of the source tables
    source_table_list = glob.glob(path.join(args.ukirt_dir,
                                            'UKIRT_year*_field*_ccd*_md.tbl'))

    # Iterating over each source table, load the set of objects detected in
    # that field pointing, CCD and year of survey
    for source_table_file in source_table_list:
        print('Scanning for known variables in source table '
              + path.basename(source_table_file))

        # Load the source table
        source_table = utils.load_ukirt_source_table(source_table_file)

        # For each catalog star, check whether it falls within the current field.
            # If it does, identify the corresponding lightcurve file, if any.
        for j,star in enumerate(var_catalog):
            if j%1000 == 0:
                print('Searching for ' + star['Name'] + ', ' + str(j)
                  + ' out of ' + str(len(var_catalog)))

            star_source_tables = find_star_in_LUT(lut, star['RA'], star['Dec'])

            if path.basename(source_table_file) in star_source_tables:
                lc_data = star['UKIRT_source_table']
                tableset = star['UKIRT_lc_files']
                if not lc_data:
                    lc_data = []
                if not tableset:
                    tableset = []

                lc = find_star_in_source_table(source_table, star['RA'], star['Dec'])

                if len(lc) > 0:
                    lc_data.append(lc)
                    tableset.append(path.basename(source_table_file))

                # If we get to this point and the lc_data list has no entries,
                # then set the entry back to None
                if len(lc_data) == 0:
                    lc_data = None
                    tableset = None

                var_catalog['UKIRT_source_table'][j] = tableset
                var_catalog['UKIRT_lc_files'][j] = lc_data


    # Output updated catalog
    utils.output_json_catalog(var_catalog, args.catalog)

def find_star_in_source_table(source_table, ra, dec):
    """
    Function to search the source table for a single field and CCD to determine
    whether at star of the given RA and Dec was detected within the field.
    If so, the function returns the year, field, CCD and star ID necessary
    to locate the corresponding lightcurve file.

    :param source_table: Path to a single source table file
    :param ra: RA of star in decimal degrees
    :param dec: Dec of star in decimal degrees
    :return: dict of the parameters needed to identify the specific lightcurve file
    """

    # Match radius tolerance:
    search_radius = (2.0 / 3600.0)*u.deg

    coord_list = SkyCoord(source_table['ra'], source_table['dec'],
                          frame='icrs', unit=(u.deg, u.deg))

    s = SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))

    sep = s.separation(coord_list)

    jdx = np.where(sep <= search_radius)[0]

    # Handle multiple matches by taking the closest
    lc_data = {}
    if len(jdx) > 0:
        lc_data = {
            'sourceid': source_table['sourceid'][jdx[0]],
            'year': source_table['obs_year'][jdx[0]],
            'bulge': source_table['bulge'][jdx[0]],
            'field': source_table['field'][jdx[0]],
            'ccdid': source_table['ccdid'][jdx[0]],
        }

    return lc_data

def find_star_in_LUT(lut, ra, dec):
    """
    Function to check whether a star is within the RA, Dec ranges of any field and CCD
    within the UKIRT survey.  If found, a list of filename(s) of the source table(s) of the
    corresponding field(s) and CCD(s) is returned.  Otherwise, an empty list is returned.

    :param ra: RA of star to search for
    :param dec: Dec of star to search for
    :return: source_tables list Filenames of source tables which may contain the star
    """

    raidx1 = np.where(ra >= lut['RA_min'])[0]
    raidx2 = np.where(ra <= lut['RA_max'])[0]
    decidx1 = np.where(dec >= lut['Dec_min'])[0]
    decidx2 = np.where(dec <= lut['Dec_max'])[0]
    idx = set(raidx1).intersection(set(raidx2))
    idx = idx.intersection(set(decidx1))
    idx = list(idx.intersection(set(decidx2)))

    return lut['source_table'][idx]

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('catalog', help='Path to input star catalog')
    parser.add_argument('lut', help='Path to UKIRT LUT')
    parser.add_argument('ukirt_dir', help='Path to directory of UKIRT source tables')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    find_ukirt_data_for_variables(args)