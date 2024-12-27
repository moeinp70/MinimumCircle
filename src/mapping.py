import pandas as pd

def preprocess_country_mapping(file_path):
    """
    Preprocess the lookup file and generate a mapping dictionary.

    Parameters:
        file_path (str): Path to the lookup file.

    Returns:
        dict: A dictionary mapping ISOCODE, UNSDCODE, and NAME0 to UNSDCODE.
    """
    file_path = "../dataset/gpw_v4_national_identifier_grid_rev11_lookup.txt"

    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep="\t")

    # Create a mapping dictionary
    country_mapping = {}
    for _, row in df.iterrows():
        isocode = row['ISOCODE']
        unsdcode = row['UNSDCODE']
        name = row['NAME0']

        if not pd.isna(isocode):
            country_mapping[isocode] = unsdcode
        if not pd.isna(unsdcode):
            country_mapping[str(unsdcode)] = unsdcode
        if not pd.isna(name):
            country_mapping[name.strip().lower()] = unsdcode



    country_mapping = {key.lower(): value for key, value in country_mapping.items()}

    print(country_mapping)

