# Software to create an index of all lightcurve files from the UKIRT
# microlensing survey of Galactic Bulge fields.

import glob
from os import path, walk
import argparse
import json

def build_index(args):
    """
    Function to build an index of all lightcurve files from the UKIRT
    microlensing survey of Galactic Bulge fields.

    The timeseries photometry for all objects within the UKIRT survey
    has been provided by the NASA Exoplanet Archive in the form of a large
    set of ASCII files, one per star.  These are organized in a series of
    directories:
    005/   006/   007/   008/   009/  010/   011/
    Within each of these is a set of subdirectories with 3-digit numbering sequence, e.g.
    000/ 001/ 002/ ... 999/
    And within each of those is a set of subdirectores with 2-digit numbers, e.g.
    00/ 01/ 02/ ... 99/
    The lightcurves for separate stars lie within these subdirectories.

    The organization pattern isn't clear here, since it doesn't seem to correspond directly to
    e.g. the fields of the survey as described in publications.  It may be organized in
    HEALpixels, arranged so that each subdirectory contains roughly the same number of
    lightcurves.

    To accelerate the retrieval of specific lightcurve files, the purpose of this
    function is to walk over the directory structure and build an index of the
    directory path to each file.

    :param args:
    :return: Output JSON catalog
    """

    # Walk over the directory tree, identifying the paths to all lightcurves
    lc_catalog = {}
    for (dirpath, subdirs, files) in walk(args.data_dir):
        if len(files) > 0:
            for f in files:
                lc_catalog[f] = path.join(dirpath,f)

    # Output the catalog
    json_object = json.dumps(lc_catalog, indent=4)
    with open(args.output_file, "w") as outfile:
        outfile.write(json_object)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Top-level directory of UKIRT lightcurve data')
    parser.add_argument('output_file', help='Path and filename of output catalog')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_index(args)