# Convenience function to extract the lightcurves from UKIRT

from os import path
import argparse
import extract_lightcurves_from_index

class Config():

    def __init__(self):
        self.var_catalog_file = None
        self.ukirt_index_dir = None
        self.year = None
        self.field = None
        self.ccd = None
        self.star_type = None
        self.output_dir = None

def run_extraction(args):
    """
    The UKIRT lightcurve archive is so large that it is most efficient to
    load the source catalog indices on a per-year, field and CCD basis,
    rather than keep re-reading these large files for each star.
    There are 296 of these index files though, so this bulk extraction
    code is designed to loop over each one in turn.
    """

    # Configuration
    years = range(2015, 2019, 1)
    fields = range(1, 65, 1)
    ccds = range(1, 5, 1)

    for yr in years:
        for field in fields:
            for ccd in ccds:
                src_table_id = 'UKIRT_year' + str(yr) + '_field' + str(field) + '_ccd' + str(ccd)
                index_file = path.join(args.ukirt_index_dir, src_table_id + '.json')
                if path.isfile(index_file):
                    print('Extracting lightcurves from ' + src_table_id)
                    conf = Config()
                    conf.var_catalog_file = args.var_catalog_file
                    conf.ukirt_index_dir = args.ukirt_index_dir
                    conf.year = str(yr)
                    conf.field = str(field)
                    conf.ccd = str(ccd)
                    conf.star_type = args.star_type
                    conf.output_dir = args.output_dir

                    extract_lightcurves_from_index.gather_data(conf)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('var_catalog_file', help='Path to variable input catalog')
    parser.add_argument('ukirt_index_dir', help='Path to top-level UKIRT lightcurve index directory')
    parser.add_argument('star_type', help='Type of variable to extract or ALL')
    parser.add_argument('output_dir', help='Path to output data directory')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    run_extraction(args)