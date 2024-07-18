# Various general-purpose utility functions

from os import path
import json
import rges_survey_definition
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column

def find_variables_in_fov(catalog, coord_type='sexigesimal'):
    """
    Function to select variables from the catalog that lie within the fields given

    :param catalog: DataFrame of known variables names, types, RA, Dec
    :return: survey_catalog: DataFrame of the variables within the survey
    """

    # Instantiate RGES survey object
    RGES = rges_survey_definition.RGESSurvey()

    # Search the catalog and find the indices of stars that lie within the survey field
    if coord_type == 'sexigesimal':
        coord_list = SkyCoord(catalog['RA'], catalog['Dec'],
                          frame='icrs', unit=(u.hourangle, u.deg))
    else:
        coord_list = SkyCoord(catalog['RA'], catalog['Dec'],
                              frame='icrs', unit=(u.deg, u.deg))

    survey_stars = RGES.find_stars_in_survey(coord_list)
    print('Identified ' + str(len(survey_stars)) + ' variables within the RGES footprint')

    # Extract the subset of the OGLE catalog for the selected stars
    survey_catalog = catalog[survey_stars]

    return survey_catalog

def output_json_catalog(catalog, output_file):
    """
    Function to output a catalog of data from an astropy Table to a JSON-format
    file.
    The first column of the catalog will be used as the dictionary key for each
    row and the remaining columns will become sub-dictionary entries.

    :param catalog: astropy table
    :param output_file: path to output JSON file
    :return: None
    """

    # Convert the table to a dictionary
    dictionary = {}
    key = catalog.colnames[0]
    columns = catalog.colnames[1:]
    for row in catalog:
        dictionary[row[key]] = {col: row[col] for col in columns}

    # Serializing to json
    json_object = json.dumps(dictionary, indent=4)

    # Output to file
    with open(output_file, "w") as outfile:
        outfile.write(json_object)

def load_json_catalog(file_path, decimal_degrees=False):
    """
    Function to load a catalog of stars in the standard format output by this code

    :param file_path:
    :param decimal_degrees: [optional, default=False] If true, convert input coordinates
                        from sexigesimal to decimal degrees
    :return: Table object of the file contents
    """

    if not path.isfile(file_path):
        raise IOError('Cannot find input catalog ' + file_path)

    with open(file_path, "r") as infile:
        json_object = json.loads(infile.read())

    names = []
    vartypes = []
    ra = []
    dec = []
    for key, entry in json_object.items():
        names.append(key)
        vartypes.append(entry['Type'])
        # Convert coordinate format if requested
        if not decimal_degrees:
            ra.append(entry['RA'])
            dec.append(entry['Dec'])
        else:
            s = SkyCoord(entry['RA'], entry['Dec'],
                         frame='icrs', unit=(u.hourangle, u.deg))
            ra.append(s.ra.deg)
            dec.append(s.dec.deg)

    # Reformat the output into a Table for ease of handling
    catalog = Table([
        Column(name='Name', data=names),
        Column(name='Type', data=vartypes),
        Column(name='RA', data=ra),
        Column(name='Dec', data=dec)
    ])

    print('Loaded ' + str(len(catalog)) + ' entries from catalog '
          + path.basename(file_path))

    return catalog