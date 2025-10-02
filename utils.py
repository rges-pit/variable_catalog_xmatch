# Various general-purpose utility functions

from os import path
import json
import rges_survey_definition
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy.io import fits
import requests
import numpy as np

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

def output_json_catalog_from_table(catalog, output_file):
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

    output_json_catalog_from_dict(dictionary, output_file)

def output_json_catalog_from_dict(catalog, output_file):
    """
    Function to output a catalog of data from an astropy Table to a JSON-format
    file.
    The first column of the catalog will be used as the dictionary key for each
    row and the remaining columns will become sub-dictionary entries.

    :param catalog: dictionary
    :param output_file: path to output JSON file
    :return: None
    """

    # Serializing to json
    json_object = json.dumps(catalog, indent=4)

    # Output to file
    with open(output_file, "w") as outfile:
        outfile.write(json_object)

    print('Output JSON catalog to ' + output_file)

def load_json_catalog_as_table(file_path, decimal_degrees=False):
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

def load_json_catalog(file_path, decimal_degrees=False):
    """
    Function to load a catalog of stars in the standard format output by this code

    :param file_path:
    :param decimal_degrees: [optional, default=False] If true, convert input coordinates
                        from sexigesimal to decimal degrees
    :return: Dictionary
    """

    if not path.isfile(file_path):
        raise IOError('Cannot find input catalog ' + file_path)

    with open(file_path, "r") as infile:
        json_object = json.loads(infile.read())

    return json_object

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

def load_ukirt_index(file_path):
    """
    Function to load the UKIRT lightcurve file index

    :param file_path:
    :return: Table of lightcurve filenames and directory paths
    """

    if not path.isfile(file_path):
        raise IOError('Cannot find UKIRT index file')

    with open(file_path, 'r') as f:
        json_object = json.loads(f.read())

    return json_object

def fetch_ogle_photometry(star_id, star_data):
    """
    Function to retrieve the OGLE timeseries photometry in multiple passbands
    for a given star in the survey's variable star catalog.  The data, if available,
    are retrieved from the survey's online archive hosted at:
    https://www.astrouw.edu.pl/ogle/ogle4/OCVS/blg/
    (note this URL is specific to the Galactic Bulge fields, suitable for RGES)

    :param star_id: Star identifier string
    :param star_data: Dictionary of star attributes
    :return:
        ilc Table or None : Timeseries photometry in I-band
        vlc Table or None : Timeseries photometry in V-band
    """

    ARCHIVE_ROOT = 'https://www.astrouw.edu.pl/ogle/ogle4/OCVS/blg/'

    # The URLS of any available timeseries photometry in each passband are
    # determinative based on the information provided by the catalog,
    # so we build these here.
    # However, the exact path to the photometry directories seems to vary
    # depending on the variability type.  Here we effectively skip variable
    # types that are unsupported.
    if star_data['Type'] == 'cephied':
        phot_dir = 'cep/phot'
    elif star_data['Type'] == 'lpv':
        phot_dir = 'lpv/phot_ogle4'
    elif star_data['Type'] == 'delta_scuti':
        phot_dir = 'dsct/phot_ogle4'
    elif star_data['Type'] == 'eclipsing_binary':
        phot_dir = 'ecl/phot_ogle4'
    elif star_data['Type'] == 'hb':
        phot_dir = 'hb/phot'
    elif star_data['Type'] == 'rrlyrae':
        phot_dir = 'rrlyr/phot'
    elif star_data['Type'] == 'cephied_type2':
        phot_dir = 't2cep/phot'
    elif star_data['Type'] == 'rot':
        phot_dir = 'rot/phot_ogle4'
    elif star_data['Type'] == 'transits':
        phot_dir = 'transits/phot_ogle4'
    elif star_data['Type'] == 'blap':
        ARCHIVE_ROOT = 'https://www.astrouw.edu.pl/ogle/ogle4/OCVS/'
        phot_dir = 'BLAP/phot/phot_ogle4'
    elif star_data['Type'] == 'cv' or star_data['Type'] == 'dn':
        ARCHIVE_ROOT = None
        phot_dir = '/data/Roman/rges_variables_wg/OGLE_lightcurves'
    else:
        phot_dir = ''

    # Access data by downloading from the OGLE online archive
    if ARCHIVE_ROOT:
        idata_url = path.join(
            ARCHIVE_ROOT,
            phot_dir, 'I',
            star_id+'.dat'
        )
        vdata_url = path.join(
            ARCHIVE_ROOT,
            phot_dir, 'V',
            star_id+'.dat'
        )

        # Check whether the corresponding lightcurve datafile is available at the
        # URLs for each bandpass.  If they are, download the data
        response = requests.get(idata_url)
        if response.status_code == 200:
            ilc = parse_ogle_lightcurve(response)
        else:
            ilc = None

        response = requests.get(vdata_url)
        if response.status_code == 200:
            vlc = parse_ogle_lightcurve(response)
        else:
            vlc = None

    # Access data from OGLE tarballs downloaded and unpacked to local disk
    # This was done for CVs/DN and the tarballs contained only I-band data
    else:
        vlc = None

        lc_file = path.join(phot_dir, star_id + '.dat')
        if path.isfile(lc_file):
            ilc = parse_ogle_lc_file(lc_file)
        else:
            ilc = None

    return ilc, vlc

def parse_ogle_lightcurve(response):
    """
    Function to parse a URL response for the lightcurve datafile of a single star

    :param response: Response object
    :return: Table of timeseries photometry
    """

    data = []
    for line in response.iter_lines():
        entries = line.decode('ascii').split()
        data.append([float(entries[0]), float(entries[1]), float(entries[2])])
    data = np.array(data)

    lc = Table(
        [
            Column(name='HJD', data=data[:,0]+2450000.0),
            Column(name='mag', data=data[:,1]),
            Column(name='mag_error', data=data[:,2])
        ]
    )

    return lc

def parse_ogle_lc_file(lc_file):
    """
    Function to parse an OGLE-format lightcurve file

    :param lc_file: str  Path to file
    :return: Table of timeseries photometry
    """

    if not path.isfile(lc_file):
        raise IOError('Cannot find OGLE lightcurve file ' + lc_file)

    data = np.loadtxt(lc_file)

    lc = Table(
        [
            Column(name='HJD', data=data[:, 0] + 2450000.0),
            Column(name='mag', data=data[:, 1]),
            Column(name='mag_error', data=data[:, 2])
        ]
    )

    return lc

def fetch_ukirt_photometry(star_id, star_data, ukirt_index, src_table_id):
    """
    Function to retrieve the UKIRT photometry for a star from a locally-mounted
    filesystem of the timeseries data from the survey.

    :param star_id: Star identifier string
    :param star_data: Dictionary of star attributes
    :param ukirt_index: Index of lightcurve files
    :param src_table_id: String identifier for working UKIRT source table
    :return: nirlc Table or None : Timeseries photometry from UKIRT
    """

    hlc = None
    klc = None

    # A star in the catalog may have zero, one or multiple UKIRT catalog entries.
    # The UKIRT sourceid can be used to look-up the file path to each
    # lightcurve from the UKIRT index.
    try:
        k = star_data['UKIRT_source_table'].index(src_table_id + '_md.tbl')
        usource = star_data['UKIRT_lc_files'][k]
        hid = usource['sourceid'] + '_h_lc.tbl'
        kid = usource['sourceid'] + '_k_lc.tbl'

        if hid in ukirt_index.keys():
            hlc = parse_ukirt_lightcurve(ukirt_index[hid])
        if kid in ukirt_index.keys():
            klc = parse_ukirt_lightcurve(ukirt_index[kid])

    except ValueError:
        pass

    return hlc, klc

def parse_ukirt_lightcurve(file_path):
    """
    Function to read and parse an UKIRT lightcurve file into an Astropy Table

    :param file_path:
    :return: Table of photometry
    """

    if path.isfile(file_path):
        with open(file_path, 'r') as f:
            file_lines = f.readlines()
            data = []
            for line in file_lines:
                if line[0:1] not in ['\\', '|']:
                    entries = line.replace('\n','').split()
                    if entries[1] not in ['null']:
                        data.append([
                            float(entries[0])+2450000.0,
                            float(entries[1]),
                            float(entries[2])
                        ])
            data = np.array(data)

            lc = Table(
                [
                    Column(name='HJD', data=data[:, 0]),
                    Column(name='mag', data=data[:, 1]),
                    Column(name='mag_error', data=data[:, 2])
                ]
            )

            return lc

    else:
        return None
def load_tess_lcs_for_star(args, tic_id, params):
    """
    Function to load one or more TESS lightcurves into a data array.
    This uses the PDCSAP_FLUX = flux series that has the common instrumental
    systematics removed.

    :param args: Arguments object
    :param tic_id: int  TESS Input Catalog identifier
    :param params: dict parameters of the star, including the filenames of any
                    lightcurves downloaded
    :return:
        lc Table Table of combined photometry from all lightcurves
        header_info: dict parameters to be added to the FITS header
    """

    if 'lc_files' not in params.keys() \
            or len(params['lc_files']) == 0:
        print('No lightcurves available to unpack for TIC ' + str(tic_id))
        return None, None

    photometry = {}
    headers = {}
    for lc_file in params['lc_files']:
        data = {'BTJD': [], 'PDCSAP_FLUX': [], 'PDCSAP_FLUX_ERR': []}
        header_info = {'sector': None, 'ccd': None}

        # Load the TESS format lightcurve
        try:
            with fits.open(path.join(args.input_dir, lc_file)) as hdul:
                header_info['sector'] = hdul[0].header['SECTOR']
                header_info['ccd'] = hdul[0].header['CCD']
                for row in hdul[1].data:
                    # Select only those lightcurve entries where the Quality flag
                    # does not indicate a spacecraft operational issue.
                    if row[9] == 0:
                        data['BTJD'].append(float(row[0]))
                        data['PDCSAP_FLUX'].append(float(row[7]))
                        data['PDCSAP_FLUX_ERR'].append(float(row[8]))

            # Convert fluxes to magnitudes.  Zeropoint derived from
            # https://heasarc.gsfc.nasa.gov/docs/tess/faq.html
            tmag = 20.44 - 2.5*np.log10(np.array(data['PDCSAP_FLUX']))
            tmag_err = (2.5 / np.log(10.0)) * np.array(data['PDCSAP_FLUX_ERR']) / np.array(data['PDCSAP_FLUX'])

            lc = Table(
                [
                    Column(name='BTJD', data=data['BTJD']),
                    Column(name='mag', data=tmag),
                    Column(name='mag_err', data=tmag_err),
                    Column(name='PDCSAP_FLUX', data=data['PDCSAP_FLUX']),
                    Column(name='PDCSAP_FLUX_ERR', data=data['PDCSAP_FLUX_ERR'])
                ]
            )
            lc.sort(['BTJD'])

            lc_idx = len(photometry)
            photometry[str(lc_idx)] = lc
            headers[str(lc_idx)] = header_info

        except OSError as e:
            print('Error loading input lightcurve ' + lc_file + ' for TIC' + str(tic_id))

    return photometry, headers

def output_multiband_lc(args, star_id, hdr, photometry):
    """
    Function to output multi-band photometry datatables as a multi-extension FITS
    binary table.

    :param star_id: Star identifier string
    :param photometry: dict of available lightcurve data Tables
    :return: None, outputs to FITS file
    """

    # Primary file header will contain basic identification information
    # for the star, plus information on the data available from the
    # lightcurve tables
    for f,lc in photometry.items():
        if lc:
            hdr['NDATA_' + f.replace('LC_','')] = len(lc)
    hdr['VSOURCE'] = 'OGLE'
    hdr['ISOURCE'] = 'OGLE'
    hdr['HSOURCE'] = 'UKIRT'
    hdr['KSOURCE'] = 'UKIRT'
    hdr['TSOURCE'] = 'TESS'

    # Add the lightcurves in each filter as a binary table extention.
    # This will create zero-length table if no data is available for a given filter.
    hdu_list = [fits.PrimaryHDU(header=hdr)]
    for f,lc in photometry.items():
        if lc:
            hdu_list.append(fits.BinTableHDU(lc, name=f))

    hdu_list = fits.HDUList(hdu_list)

    # Output to disk:
    file_path = get_lc_path(args, star_id)
    if not path.isdir(path.dirname(file_path)):
        makedirs(path.dirname(file_path))
    hdu_list.writeto(file_path, overwrite=True)

    return file_path

def make_lc_header(star_id, star_data):
    """
    Function to generate a standardized lightcurve file header
    :param star_id:
    :param star_data:
    :return: FITS header object
    """

    hdr = fits.Header()
    hdr['NAME'] = star_id
    hdr['RA'] = star_data['RA']
    hdr['DEC'] = star_data['Dec']
    hdr['VARTYPE'] = star_data['Type']

    return hdr


def get_lc_path(args, star_id):
    """
    Function to return a standardized path to a lightcurve file
    :param args: Program commandline arguments
    :param star_id: Identifier for star
    :return:
    :param file_path: string, path to lightcurve
    """

    return path.join(args.output_dir, star_id + '_multiband_lc.fits')

def load_multiband_lc(file_path):
    """
    Function to read a multiband lightcurve file

    :param file_path: String path to input file
    :return:
    :param header: Dictionary of lightcurve header parameters
    :param lightcurves: Dictionary of lightcurves, as astropy Tables
    """

    if not path.isfile(file_path):
        raise IOError('Cannot find lightcurve file at ' + file_path)

    with fits.open(file_path) as hdu_list:
        header = hdu_list[0].header
        lightcurves = {}
        for hdu in hdu_list:
            if 'LC_' in hdu.name:
                data = [[dp[0], dp[1], dp[2]] for dp in hdu.data]
                data = np.array(data)

                lc = Table([
                    Column(name='HJD', data=data[:, 0]),
                    Column(name='mag', data=data[:, 1]),
                    Column(name='mag_error', data=data[:, 2])
                ])
                lightcurves[hdu.name] = lc

    return header, lightcurves