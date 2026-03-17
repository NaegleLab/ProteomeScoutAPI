import os
import json

#grab package directory (include specific folder)
API_DIR = os.path.dirname(os.path.abspath(__file__))
FIGSHARE_ID = 31129045
###### Customizable configuration ####################
CONFIG_FILE = API_DIR + "/configuration.json"
with open(CONFIG_FILE, 'r') as f:
    configuration = json.load(f)


VERSION = configuration["DATASET_VERSION"]
UPDATE = configuration["UPDATE"]
DATASET_DIR = configuration["DATASET_DIR"]
if DATASET_DIR is None:
    DATASET_DIR = API_DIR


def update_configuration(dataset_dir = None, version = None, update = None):
    """Update configuration file for ProteomeScoutAPI.

    Parameters
    ----------
    dataset_dir : str, optional
        Path to the dataset directory.
    version : int, optional
        Version number of the dataset.
    update : bool, optional
        Whether to force an update of the dataset if available or if different version is requested.
    """
    global VERSION, DATASET_DIR, UPDATE
    if dataset_dir is not None:
        DATASET_DIR = dataset_dir
    if version is not None:
        VERSION = version
    if update is not None:
        UPDATE = update

    config_data = {
        "DATASET_VERSION": VERSION,
        "DATASET_DIR": DATASET_DIR,
        "UPDATE": UPDATE
    }
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config_data, f, indent=4)




