# Various general-purpose utility functions
import json

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
    print(dictionary)

    # Serializing to json
    json_object = json.dumps(dictionary, indent=4)

    # Output to file
    with open(output_file, "w") as outfile:
        outfile.write(json_object)