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

    print('Output JSON catalog to ' + output_file)

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

def load_ukirt_source_table(file_path):
    """
    Function to load a source table file from the UKIRT microlensing survey.
    Each source table provides the metadata for stars detected within
    a single field pointing and one of the four CCDs.

    Based on the header information from the source tables, the columns included are:
    sourceid        Mission specific source identification
    phot_method     Data pipeline used to produce the photometry
    obs_year        Calendar year of observations
    bulge           Bulge
    field           Time Series Minimum integer time
    ccdid           Time Series Maximum integer time
    hjdstart (days) Time Series Minimum Julian Date
    hjdstop (days)  Time Series Maximum Julian Date
    hjd_ref (days)  Base Julian Date
    ra (degrees)    Object Right Ascension
    dec (degrees)   Object Declination
    h_mag (mag)     H Magnitude
    j_mag (mag)     J Magnitude
    k_mag (mag)     K Magnitude
    npts            Points in Light Curve
    k2c9_flag       Flag indicating overlap with K2 C9 field
    ukirt_evt_flag  Is this a ukirt event
    ukirt_id        UKIRT survey event identification string
    ogle_evt_flag   Is this an ogle event
    ogle_id         OGLE survey event identification string
    moa_evt_flag    Is this a MOA event
    moa_id          MOA survey event identification string
    moa_star_id     MOA survey event label string
    statnpts        number of points used in MAG statistics calculations
    minvalue (mag)  minimum value of MAG column
    maxvalue (mag)  maximum value of MAG column
    mean (mag)      mean value of MAG column
    stddevwrtmean (mag) Standard deviation with respect to mean of MAG column
    median (mag)    median value of MAG column
    stddevwrtmed (mag)  Standard deviation with respect to median of MAG column
    n5sigma         Number of points beyond 5 stddev wrt MAG median
    f5sigma         Fraction of points beyond 5 stddev wrt MAG median
    medabsdev (mag) Median absolute deviation of MAG column
    chisquared      Reduced Chi Squared wrt MAG median
    range595 (mag)  Range of MAG, excluding 5% of minimum and 5% of maximum

    :param file_path: string Path to file
    :return: Table of file contents
    """

    COLUMN_LIST = [
        'sourceid',
        'phot_method',
        'obs_year',
        'bulge',
        'field',
        'ccdid',
        'hjdstart',
        'hjdstop',
        'hjd_ref',
        'ra',
        'dec',
        'h_mag',
        'j_mag',
        'k_mag',
        'npts',
        'k2c9_flag',
        'ukirt_evt_flag',
        'ukirt_id',
        'ogle_evt_flag',
        'ogle_id',
        'moa_evt_flag',
        'moa_id',
        'moa_star_id',
        'statnpts',
        'minvalue',
        'maxvalue',
        'mean',
        'stddevwrtmean',
        'median',
        'stddevwrtmed',
        'n5sigma',
        'f5sigma',
        'medabsdev',
        'chisquared',
        'range595',
    ]
    FLOAT_COLUMNS = ['ra', 'dec']

    if not path.isfile(file_path):
        raise IOError('Cannot find UKIRT source table ' + file_path)

    data = {key: [] for key in COLUMN_LIST}

    with open(file_path, 'r') as f:
        file_lines = f.readlines()

        for line in file_lines:
            if line[0:1] not in ["\\", '|']:
                entries = line.replace('\n','').split()
                if len(entries) == len(COLUMN_LIST):
                    [
                        data[key].append(float(entries[i]) if key in FLOAT_COLUMNS else entries[i])
                        for i,key in enumerate(COLUMN_LIST)
                    ]

    # Reformat the output into a Table for ease of handling
    table_columns = [Column(name=key, data=data[key]) for key in COLUMN_LIST]
    source_table = Table(table_columns)

    return source_table

def load_ukirt_lut(lut_file):
    """
    Function to load the UKIRT Look-Up Table

    :param lut_file:
    :return: Table of LUT contents
    """

    if not path.isfile(lut_file):
        raise IOError('Cannot find UKIRT survey LUT file')

    with open(lut_file, 'r') as f:
        json_object = json.loads(f.read())

    source_tables = []
    ramin = []
    ramax = []
    decmin = []
    decmax = []
    for key, entry in json_object.items():
        source_tables.append(key)
        ramin.append(float(entry['ra_min']))
        ramax.append(float(entry['ra_max']))
        decmin.append(float(entry['dec_min']))
        decmax.append(float(entry['dec_max']))

    # Reformat the output into a Table for ease of handling
    lut = Table([
        Column(name='source_table', data=source_tables),
        Column(name='RA_min', data=ramin),
        Column(name='RA_max', data=ramax),
        Column(name='Dec_min', data=decmin),
        Column(name='Dec_max', data=decmax)
    ])

    print('Loaded ' + str(len(lut)) + ' entries from catalog '
          + path.basename(lut_file))

    return lut