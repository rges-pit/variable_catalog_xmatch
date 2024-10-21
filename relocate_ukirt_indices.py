# Code to update the file paths to the table files of the UKIRT catalog index files

from os import path, rename
import glob
import argparse
import utils

# Fetch commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('index_dir', help='Path to directory containing UKIRT index files')
parser.add_argument('old_prefix', help='Path prefix to old location of the UKIRT source catalog files')
parser.add_argument('new_prefix', help='Path prefix to new location of the UKIRT source catalog files')
args = parser.parse_args()

# Loop over all JSON-format catalog index files, updating the route file path given for
# table file
index_files = glob.glob(path.join(args.index_dir, '*.json'))

for file_path in index_files:
    index = utils.load_ukirt_index(file_path)

    # Update the catalog table file paths in the index to the path on local disk
    # e.g. /Volumes/U/UKIRT/ -> /media/rstreet/U/UKIRT/
    new_index = {}
    for key, cat_path in index.items():
        new_index[key] = cat_path.replace(args.old_prefix, args.new_prefix)

    old_file = file_path.replace('.json', '_old.json')
    rename(file_path, old_file)

    utils.output_json_catalog_from_dict(new_index, file_path)

    print('Relocated entries for ' + file_path)