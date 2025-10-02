from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
import argparse
from os import path
import utils
import identify_known_ogle_variables

def add_var_type(args):
    """
    Function to add information about a new type of variables to an existing catalog of variable stars.
    """

    # Load the configuration info on the variable type selected.
    var_configs = identify_known_ogle_variables.load_config()
    if args.var_type not in var_configs['variable_catalogs'].keys():
        raise IOError(
            'No configuration for variable type ' + args.var_type
            + ' available in identify_known_ogle_variables.load_config'
          )
    config = var_configs['variable_catalogs'][args.var_type]

    # Load the existing variable star catalog
    var_catalog = utils.load_json_catalog(args.var_catalog_file)

    # Load the list of new variables
    new_catalog = {
        'name': [],
        'type': [],
        'ra': [],
        'dec': []
    }
    new_catalog = identify_known_ogle_variables.parse_ogle_variable_catalog(
        new_catalog,
        var_configs['data_dir'],
        args.var_type,
        config
    )
    ogle_catalog = Table([
        Column(name='Name', data=new_catalog['name']),
        Column(name='Type', data=new_catalog['type']),
        Column(name='RA', data=new_catalog['ra']),
        Column(name='Dec', data=new_catalog['dec'])
    ])

    # Identify those objects within the field of view of all given fields
    survey_catalog = utils.find_variables_in_fov(ogle_catalog, coord_type='sexigesimal')

    # Add the new variables to the var_catalog
    for star in survey_catalog:
        if star['Name'] not in var_catalog.keys():
            coord = SkyCoord(star['RA'], star['Dec'], frame='icrs', unit=(u.hourangle, u.deg))
            var_catalog[star['Name']] =  {
                "Type": star['Type'],
                "RA": coord.ra.deg,
                "Dec": coord.dec.deg,
                "VVV_ID": None,
                "VVV_type": None,
                "Gaia_ID": None,
                "Gaia_type": None,
                "UKIRT_source_table": None,
                "UKIRT_lc_files": None
            }

    # Output updated version of the variable catalog
    utils.output_json_catalog_from_dict(var_catalog, args.var_catalog_file)

def get_args():
    """
        Gather user-provided input

        :return: args dict Set of argument key-value pairs
        """

    parser = argparse.ArgumentParser()
    parser.add_argument('var_catalog_file', help='Path to variable catalog file')
    parser.add_argument('var_type', help='Type of new variable to add')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    add_var_type(args)