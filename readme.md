# **ProteomeScoutAPI**


A lightweight API to talk to ProteomeScout flatfiles in Python and annotate phosphoproteomic datasets with protein and site-specific information

See full [documentation](https://naeglelab.github.io/ProteomeScoutAPI/) for installation and use cases.

## About ProteomeScout
Version 3.0 - March 2026

Current Version: Naegle Lab, University of Virginia [https://github.com/naegleLab](https://github.com/naegleLab)

Originally: By Alex Holehouse, Washington University in St. Louis.

## Overview
ProteomeScoutAPI is a Python module which can be used to connect to and parse ProteomeScout flatfiles. Specifically, the goal of this module is to allow anyone to interact with ProteomeScout data without the need to

1. Repeatedly query the ProteomeScout sever

2. Have any knowledge of SQL, or use an SQL-Python ORM

3. Facilitate rapid exploration of the ProteomeScout dataset



## Available Information

The ProteomeScout flat file contains information associated with all proteins in the human proteome, including:
- Protein sequence
- Protein domains (from InterPro or UniProt)
- Protein structures (alpha helices, beta strands, etc)
- Macromolecular structures (disorder, etc.)
- PTMs
- Gene Ontology terms
- Exon information
- Kinase activation loops
- SpY-C predictions (likelihood a sequence is recognized by an SH2 domain)

## Usage

ProteomeScoutAPI allows for fast querying of protein-specific information, allowing for users to:
1. Explore properties of individual protein(s) of interest (domains, PTMs, structure, etc.)
2. Annotate a phosphoproteomic dataset with additional information about the proteins and sites in the dataset (e.g. domains, GO terms, etc.)
       
 


# Contributing code
The code here is incredibly simple, and the few methods presented give a good example of how one should parse the ProteomeScout records. 

For step-by-step guidance on adding a new field so it propagates through API annotations, species datasets, and proteomic dataset annotation, see [FEATURE_INTEGRATION_README.md](FEATURE_INTEGRATION_README.md).

Trying to use your own version of the database file? For example, testing a new feature integration. If so, you will want to make sure to force the API to point to and NOT update to the latest Figshare version. 

# Dataset root now that you added ProteomeScout_Dataset under April_cleaned.
``` dataset_root = Path("<your/local/path>")  # <-- UPDATE THIS TO YOUR LOCAL PATH
expected_data_file = dataset_root / "ProteomeScout_Dataset" / "data.tsv"
if not expected_data_file.exists():
    raise FileNotFoundError(f"Could not find local dataset file: {expected_data_file}")
config.UPDATE = False  # Keep this False so version mismatch checks do not auto-download.
config.DATASET_DIR = str(dataset_root)
API = ProteomeScoutAPI(version=set_to_your_version_indicated_in_json, update=False)  # Keep this
```