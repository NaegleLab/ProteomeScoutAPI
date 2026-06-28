# ProteomeScoutAPI
#
# ===============================================================================
# ABOUT
# ===============================================================================
# Version 2.0
# 
# June 2020
#
# Oritinally By Alex Holehouse, Washington University in St. Louis 
# Contact alex.holehouse@gmail.com 
# This version maintained by Naegle Lab
# Contact Kristen Naegle at kmn4mj@virginia.edu
#
#
#
# ===============================================================================
# OVERVIEW
# ===============================================================================
# ProteomeScoutAPI is a Python module which can be used to connect to and parse
# ProteomeScout flatfiles. 
#
# Specifically, the goal of this module is to allow anyone to interact with 
# ProteomeScout data without the need to
#
# 1) Repeatedly query the ProteomeScout sever
#
# 2) Have any knowledge of SQL, or use an SQL-Python ORM
#
# 3) Facilitate rapid exploration of the ProteomeScout dataset
#
# ===============================================================================
# Usage
# ===============================================================================
# The general approach to usage is as follows
#
# 0) Import the proteomeScoutAPI module
#  
# 1) Create the API object by loading a flat file. ProteomeScout flat files can 
#    be obtained from https://proteomescout.wustl.edu/compendia
#
#    PTM_API = ProteomeScoutAPI('<flat filename goes here>')
#
# 2) Query the API object using the functions. Available functions are
#
#    PTM_API.get_phosphosites(ID)
#    PTM_API.get_PTMs(ID)
#    PTM_API.get_mutations(ID)
#
#    For more information on these functions I suggest reading the rest of the
#    source code, or once you've loaded an API object type
#
#    help(PTM_API.get_mutations)
#
# 3) The PTM_APU.uniqueKeys is a list of the unique accession numbers to provide
#    an easy way to loop over all the unique entries. NOTE THAT IDS HAVE REDUNDANCY
#    which is deliberate (and makes interfacing easy) but cycling over all the IDs
#    in the API object would be incorrect and lead to double counting
#
# ===============================================================================
# EXAMPLES
# ===============================================================================
# from proteomeScoutAPI import ProteomeScoutAPI 
#
# 
# ID = "O46631"
#  
# PTM_API = ProteomeScoutAPI("proteomescout_mammalia_20140831.tsv")
#
# PTM_API.get_mutations(ID)
#
# PTM_API.get_PTMs(ID)
#
# # the following loop prints all the phosphosites in all the proteins
# for ID in PTM_API.uniqueKeys():
#    print PTM_API.get_phosphosites(ID)
#
# ===============================================================================
import os, json, contextlib, re
import pandas as pd
import numpy as np


# import helprs
from proteomeScoutAPI import figshare, config, helpers



# Exception if file is bad
class BadProteomeScoutFile(Exception):
    pass

class ProteomeScoutAPI:

    """
    Class for interacting with ProteomeScout flat files.

    Parameters
    ----------
    version : int, optional
        Version number of the dataset to use. Defaults to the version specified in the configuration. If None, uses the latest version.
    update : bool, optional
        Whether to check for updates to the dataset upon initialization. Defaults to the value specified in the configuration (True unless changed).
    
    Attributes
    ----------
    database : dict
        Internal database storing ProteomeScout data.
    uniqueKeys : list
        List of unique accession numbers in the dataset.
    version : int
        The version number of the dataset being used.
    
    """

    POSITION_ANNOTATION_REGISTRY = [
        {
            'context_key': 'interpro',
            'name_column': 'Domain_Names_InterPro',
            'extra_column': 'InterPro_IDs',
            'include_boolean': False,
        },
        {
            'context_key': 'uniprot',
            'name_column': 'Domain_Names_UniProt',
            'include_boolean': False,
        },
        {
            'context_key': 'structure',
            'name_column': 'Structures',
            'include_boolean': False,
        },
        {
            'context_key': 'macro',
            'name_column': 'Macro_Molecular_Structures',
            'include_boolean': False,
        },
        {
            'context_key': 'activation_loop',
            'name_column': 'Activation_Loop_Quality',
            'boolean_column': 'In_Activation_Loop',
            'include_boolean': True,
        },
        {
            'context_key': 'exons',
            'name_column': 'Exons',
            'extra_column': 'Exon_Constitutive',
            'boolean_column': 'In_Exon',
            'include_boolean': True,
        },
    ]
    def __init__(self, version = config.VERSION, update = config.UPDATE):
        #initialize figshare interface       
        try:
            self.figshare_interface = figshare.FigshareInterface(config.FIGSHARE_ID)
            #get latest version
            self.latest_version = self.figshare_interface.get_latest_dataset_version()
            #set version number, default to latest if None
            self.version = version if version is not None else self.latest_version
        except Exception as e:
            print(f"Unable to connect to FigShare to assess available versions (Error: {e}). Will rely on local version if available.")
            self.version = version
            self.figshare_interface = None
            self.latest_version = None

        #look in package directory for dataset
        self.dataset_dir = config.DATASET_DIR
        self.update = update

        
        #check if data.tsv file exists in dataset_dir
        self.datafile = os.path.join(self.dataset_dir, "ProteomeScout_Dataset", "data.tsv")
        self.meta_file = os.path.join(self.dataset_dir, "ProteomeScout_Dataset", "metadata.json")
        if not os.path.isfile(self.datafile):
            if self.figshare_interface is not None:
                print("ProteomeScout dataset file not found. Downloading data to %s"%(self.dataset_dir))
                self.download_data(version = self.version)
            else:
                raise RuntimeError("No ProteomeScout data file found locally, and unable to connect to FigShare to download. Please try again later.")
        else:
            self.__checkDatasetVersion(update = self.update)

    

        self.database={}
        self.uniqueKeys=[]
        # this will throw and exception if there's a problem with the file
        self.__checkFile(self.datafile)
        
        self.__buildAPI(self.datafile)

        #load citations
        citations_file = os.path.join(self.dataset_dir, "ProteomeScout_Dataset", "citations.tsv")
        self.citations = pd.read_csv(citations_file, sep='\t')
        self.experiment_current_flags = self.__build_experiment_current_flags(self.citations)

    def __build_experiment_current_flags(self, citations_df):
        """
        Build a fast lookup for experiment visibility keyed by experiment ID.

        Parameters
        ----------
        citations_df : pd.DataFrame
            Citation metadata loaded from citations.tsv.

        Returns
        -------
        dict
            Mapping of experiment ID string to whether the experiment is current.
        """
        experiment_flags = {}

        if citations_df is None or citations_df.empty:
            return experiment_flags

        if 'Experiment ID' not in citations_df.columns or 'Current' not in citations_df.columns:
            return experiment_flags

        for _, row in citations_df[['Experiment ID', 'Current']].dropna(subset=['Experiment ID']).iterrows():
            experiment_id = str(row['Experiment ID']).strip()
            if experiment_id.endswith('.0'):
                experiment_id = experiment_id[:-2]

            current_value = str(row['Current']).strip().lower()
            experiment_flags[experiment_id] = current_value in {'1', 'true', 'yes'}

        return experiment_flags

    def __parse_evidence_tokens(self, evidence_token):
        """
        Extract experiment IDs from one PTM evidence token.

        Parameters
        ----------
        evidence_token : str
            One semicolon-delimited evidence entry.

        Returns
        -------
        list of str
            Experiment IDs found in the token.
        """
        if evidence_token is None:
            return []

        return re.findall(r'\d+', str(evidence_token))

    def __filter_modifications_by_visibility(self, record, include_hidden=False):
        """
        Filter PTMs so hidden-only evidence is excluded by default.

        Parameters
        ----------
        record : dict
            Protein record from the in-memory database.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.

        Returns
        -------
        list of tuple
            Filtered PTM tuples as (position, residue, modification-type).
        list of str
            Evidence tokens aligned to the filtered PTMs.
        """
        mods = record.get("modifications", "")
        if not mods or mods.strip() == "":
            return [], []

        mods_clean = helpers.clean_PTM_string(mods)
        evidence = record.get("evidence", "")
        evidence_tokens = [token.strip() for token in str(evidence).split(';')]

        if include_hidden or not evidence_tokens or len(evidence_tokens) != len(mods_clean):
            trimmed_evidence = evidence_tokens[:len(mods_clean)]
            if len(trimmed_evidence) < len(mods_clean):
                trimmed_evidence.extend([''] * (len(mods_clean) - len(trimmed_evidence)))
            return mods_clean, trimmed_evidence

        filtered_mods = []
        filtered_evidence = []
        for mod, evidence_token in zip(mods_clean, evidence_tokens):
            experiment_ids = self.__parse_evidence_tokens(evidence_token)
            if not experiment_ids:
                filtered_mods.append(mod)
                filtered_evidence.append(evidence_token)
                continue

            flags = [self.experiment_current_flags.get(experiment_id) for experiment_id in experiment_ids]
            known_flags = [flag for flag in flags if flag is not None]

            if not known_flags or any(known_flags):
                filtered_mods.append(mod)
                filtered_evidence.append(evidence_token)

        return filtered_mods, filtered_evidence

    def __checkDatasetVersion(self,update = False):
        """
        During initialization, check if the dataset version found in the dataset directory matches the requested version.
        
        Parameters
        ----------
        update : bool, optional
            Whether to update the dataset if the versions do not match. Defaults to False.
        """
        #check to make sure the dataset is the expected version
        if self.version is None:
            print('No version specified and could not assess available versions. Proceeding without version check.')
            return
        elif os.path.isfile(self.meta_file):
            #load metadata to get version
            with open(self.meta_file, 'r') as f:
                self.metadata = json.load(f)

            #check 
            current_version = self.metadata.get('version_number', None)
            #if version is indicated, check if it matches the requested version
            if self.version is not None:
                if current_version != self.version:
                    print(f"WARNING: Requested version {self.version} does not match existing dataset version {current_version}.")
                    if update:
                        print("Updating dataset to requested version...")
                        self.download_data(version = self.version)
                    else:
                        print("Proceeding with existing dataset version. To update to the latest version, use download_data() or set the update parameter to True.")
                        self.version = current_version
            else:
                self.check_for_updates(update = update)
        else:
            print("WARNING: Data file found, but no metadata file found. Proceeding without version check.")

    def __checkFile(self, filename):
        """
        Internal function to check if there is have a valid proteomescout flat file.

        Parameters
        ----------
        filename : str
            Path to the ProteomeScout flat file to check.

        Raises
        ------
        BadProteomeScoutFile
            If the file is not a valid ProteomeScout flat file.
        """
        
        try:
            with open(filename, 'r') as f:
                first_line = f.readline().strip()

            headers = [header.strip() for header in first_line.split("\t") if header.strip()]
            required_headers = {
                "protein_id", "accessions", "acc_gene", "protein_name", "species",
                "sequence", "modifications", "evidence", "uniprot_domains",
                "macro_molecular", "structure", "GO_terms", "uniprot_id",
                "updated", "error_code", "Interpro_domains", "swissprot_nr"
            }

            if len(headers) < 17 or not required_headers.issubset(set(headers)):
                raise BadProteomeScoutFile("N/A")
            
            optional_headers = {"activation_loop", "exons", "spyc_predictions"}
            missing_optional = optional_headers - set(headers)
            if missing_optional:
                print(f"WARNING: The following optional headers are missing from the API: {', '.join(missing_optional)}. This likely means you are using an older version of the API dataset. If you would like to use these, please update the dataset with `update_to_latest()` and/or set the update parameter to True.")
            
                
        except:
            raise BadProteomeScoutFile("Invalid ProteomeScout flat file %s.\nFile is invalid or corrupted" % str(filename))
        

    def __buildAPI(self, datafile):
        """
        Internal function to build the API database from the ProteomeScout flat file.
        
        Parameters
        ----------
        datafile : str
            Path to the ProteomeScout flat file to parse.

        Postconditions
        --------------
        Populates self.database and self.uniqueKeys with data from the flat file.
        """
        # read in the file line by line
        with open(datafile) as f:
            content = f.readlines()
            

        
        # read the flatfile headers and parse
        # the contents
        headers_raw = content[0].strip().split('\t')

        headers=[]
        for i in headers_raw:
            headers.append(i.strip())

        try:
            uniprot_id_index = headers.index("uniprot_id")
        except ValueError:
            uniprot_id_index = 12

        # remove the header line from the read in data
        content.pop(0)

        
        # for each line in the data
        for line in content:

            # split the record and get the canonical ID
            record          = line.strip().split('\t')
            
            # Skip empty lines
            if not record or not record[0].strip():
                continue
            
            ProteomeScoutID = record[0]
            
            # Use uniprot_id as the primary lookup key
            uniprot_id = record[uniprot_id_index].strip() if len(record) > uniprot_id_index and record[uniprot_id_index].strip() else ProteomeScoutID
                
            # add to the list of unique keys
            self.uniqueKeys.append(uniprot_id)

            # now construct the object dictionary
            OBJ={}
            for i in range(1, len(headers)):
                OBJ[headers[i]] = record[i].strip() if i < len(record) else ""

            # ALWAYS add the record via uniprot_id (or ProteomeScoutID if no uniprot_id)
            if uniprot_id in self.database:
                # if the newly found record has more PTMs associated with it than the 
                # original record then overwrite 
                if len(OBJ['modifications'].split(";")) > len(self.database[uniprot_id]['modifications'].split(";")):
                    self.database[uniprot_id] = OBJ
            else:
                self.database[uniprot_id] = OBJ

    def check_for_updates(self, update=False):
        """
        Check if a newer version of the ProteomeScout dataset is available on FigShare.

        Parameters
        ----------
        update : bool, optional
            Whether to update the dataset if a newer version is available. Defaults to False.
        
        Postconditions
        --------------
        If a newer version is available, prints statement indicating the newest version available. If update is True, downloads and updates to the latest version.
        """
        if self.version is None and self.latest_version is None:
            print('Unable to check for updates because unable to connect to FigShare.')
        elif self.version is None or self.latest_version > self.version:
            print(f"A newer version ({self.latest_version}) of the ProteomeScout dataset is available. You are using version {self.version}. Run update_to_latest() to download the latest version.")
            if update:
                self.__update()
        else:
            print("You are using the latest version of the ProteomeScout dataset.")
            return None
        
    def update_to_latest(self):
        """
        Download and update to the latest version of the ProteomeScout dataset from FigShare.

        Postconditions
        --------------
        If a newer version is available, downloads and updates the dataset to the latest version, rebuilding the API with the new dataset.
        """
  
        if self.latest_version > self.version:
            print(f"Updating to version {self.latest_version} of the ProteomeScout dataset...")
            self.download_data(version = self.latest_version)
            self.version = self.latest_version
            # Rebuild the API with the new dataset
            self.database={}
            self.uniqueKeys=[]
            self.__checkFile(self.datafile)
            self.__buildAPI(self.datafile)
            print("Update complete.")
        else:
            print("You are already using the latest version of the ProteomeScout dataset.")



    def download_data(self, version = None):
        """
        Retrieves proteomescout data files that are the companion for this version release from FigShare. Will download and decompress the files to self.dataset_dir in a "ProteomeScout_Dataset" folder

        Parameters
        ----------
        version : int, optional
            Version number of the dataset to download from FigShare. If None, uses the latest version.

        Postconditions
        --------------
        Downloads and decompresses the ProteomeScout dataset to self.dataset_dir/ProteomeScout_Dataset        
        """
        #download files
        if version is None:
            version = self.version
        self.figshare_interface.download(version, self.dataset_dir)

        #save version number
        self.version = version

    def get_PTMs(self, ID, output_format='list', include_hidden=False):
        """
        Return all PTMs associated with the ID in question.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        list of tuples, pd.DataFrame, or int
            tuple for each PTM associated with the ID(position, residue, modification-type). If output_format is 'table', returns a pandas DataFrame with columns ['Position', 'Residue', 'Modification_Type'].
            Returns an empty list if no modifications are found. 
        
        Postconditions
        --------------
        If ID is not found, returns -1
        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        mods_clean, _ = self.__filter_modifications_by_visibility(record, include_hidden=include_hidden)
        if not mods_clean:
            if output_format == 'table':
                return pd.DataFrame(columns=['Position', 'Residue', 'Modification_Type'])
            return []

        if output_format == 'table':
            # Convert to DataFrame
            mods_clean = pd.DataFrame(mods_clean, columns=['Position', 'Residue', 'Modification_Type'])
        return mods_clean
    
    def get_structure(self, ID, output_format = 'list'):
        """
        Return all structures associated with the ID in question.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.
        
        Returns
        -------
        list of tuples, pd.DataFrame, or int
            tuple for each structures associated with the ID (domain_name, start_position, end_position). If output_format is 'table', returns a pandas DataFrame with columns ['Structure_Name', 'Start_Position', 'End_Position'].
            Returns returns an empty list or dataframe if no structures are found. 
       

        Postconditions
        --------------
        If ID is not found or there is a problem with the start and end position, returns -1  
        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")
        try:
            record = self.database[ID]
        except KeyError:
            return -1
            
        structs = record["structure"]
        
        if not structs or structs.strip() == "":
            return []
               
        structs_raw=structs.split(";")
        structs_clean =[]
        for i in structs_raw:
            if i:
                tmp = i.strip()
                tmp = tmp.split(":")
                if len(tmp) >= 2:
                    structs_clean.append((tmp[0], tmp[1], tmp[2]))

        if output_format == 'table':
            # Convert to DataFrame
            structs_clean = pd.DataFrame(structs_clean, columns=['Structure_Name', 'Start_Position', 'End_Position'])
        return structs_clean
    
    def get_macro_molecular(self, ID, output_format = 'list'):
        """
        Return all macro-molecular structures associated with the ID in question.

        Returns
        -------
        list of tuples, pd.DataFrame, or int
            tuple for each macro-molecular structure associated with the ID (structure_name, start_position, end_position). If output_format is 'table', returns a pandas DataFrame with columns ['Macro_Name', 'Start_Position', 'End_Position'].
            Returns returns an empty list if no macro-molecular structures are found.

        POSTCONDITIONS
        --------------
        If ID is not found, returns -1
        
        Removes any malformed entries that don't have name:start:stop


        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
            
        macro_mol = record["macro_molecular"]
        
        if not macro_mol or macro_mol.strip() == "":
            return []
               
        macro_mol_raw = macro_mol.split(";")
        macro_mol_clean = []
        for i in macro_mol_raw:
            if i:
                tmp = i.strip()
                parts = tmp.split(":")
                if len(parts) >= 3:
                    name = parts[0]
                    start = parts[1]
                    stop = parts[2]
                    try:
                        int(start)
                        int(stop)
                    except ValueError:
                        continue
                    macro_mol_clean.append((name, start, stop))
                else:
                    print("ERROR: the macro-molecular structure did not match expected format %s"%(i))

        if output_format == 'table':
            macro_mol_clean = pd.DataFrame(macro_mol_clean, columns=['Macro_Name', 'Start_Position', 'End_Position'])
        return macro_mol_clean
    
    def get_domains(self, ID, domain_type=None, output_format = 'list'):
        """
        Return all domains associated with the ID in question.
        For interpro domains domain_type is 'interpro'
        For UniProt domains domain_type is 'uniprot'
        If domain_type is not specified, returns a dictionary with both.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        domain_type : str, optional
            Type of domain to retrieve ('interpro' or 'uniprot'). If None, returns both types.
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.

        Returns
        -------
        list of tuples, dict, int, or pandas.DataFrame
            tuples for each domain associated with the ID (domain_name, start_position, end_position, domain_id). Interpro domains include the interpro_id as the last element in the tuple; uniprot domains have None as the last element. If output_format is 'table', returns a pandas DataFrame with columns ['Domain_Name', 'Start_Position', 'End_Position', 'Domain_ID'] instead of list.

            Returns an empty list if no domains are found

            Returns a dictionary with 'interpro' and 'uniprot' keys if domain_type is not specified.

        Postconditions
        --------------
        Returns -1 if unable to find the ID

        Returns -2 if an unrecognized domain_type is specified

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        # If no domain_type specified, return both
        if domain_type is None:
            return {
                'interpro': self.get_domains(ID, 'interpro', output_format=output_format),
                'uniprot': self.get_domains(ID, 'uniprot', output_format=output_format)
            }
        
        domain_type = domain_type.lower()
        if domain_type == 'interpro':
            doms = record["Interpro_domains"]
        elif domain_type == 'uniprot':
            doms = record["uniprot_domains"]
        else:
            print("%s is an unrecognized domain type. Use 'interpro' or 'uniprot'"%(domain_type))
            return -2

        #check for atypical records and replace the extra ; to avoid errors in parsing
        doms = doms.replace("; atypical", " atypical")
        doms = doms.replace("; truncated", " truncated")
        doms = doms.replace("; second part", " second part")
        doms = doms.replace("; first part", " first part")
        
        doms_raw=doms.split(";")
        doms_clean =[]
        for i in doms_raw:
            if i:
                tmp = i.strip()
                if domain_type == 'interpro':
                    # interpro format: name:interpro_id:start:stop
                    parts = tmp.split(":")
                    if len(parts) >= 4:
                        name = parts[0]
                        interpro_id = parts[1]
                        start = parts[2]
                        stop = parts[3]
                        doms_clean.append((name, start, stop, interpro_id))
                    else:
                        print("ERROR: the interpro domain did not match expected format %s"%(i))
                else:
                    # uniprot format: name:start:stop
                    parts = tmp.split(":")
                    if len(parts) >= 3:
                        name = parts[0]
                        start = parts[1]
                        stop = parts[2]
                        doms_clean.append((name, start, stop, None))
                    else:
                        print("ERROR: the uniprot domain did not match expected format %s"%(i))

        if output_format == 'table':
            # Convert to DataFrame
            doms_clean = pd.DataFrame(doms_clean, columns=['Domain_Name', 'Start_Position', 'End_Position', 'Domain_ID'])

        return doms_clean

    def get_activation_loops(self, ID, output_format='list'):
        """
        Return activation loop regions associated with the ID.

        Activation loop entries are expected in the format start:stop:quality,
        separated by semicolons.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.

        Returns
        -------
        list of tuples, pd.DataFrame, or int
            List of tuples in the form (quality, start_position, end_position).
            If output_format is 'table', returns a DataFrame with columns
            ['Quality', 'Start_Position', 'End_Position'].

            Returns an empty list when activation_loop data are absent.

        Postconditions
        --------------
        If ID is not found, returns -1
        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")

        try:
            record = self.database[ID]
        except KeyError:
            return -1

        loops = record.get("activation_loop", "")
        if not loops or loops.strip() == "":
            return []

        loops_raw = loops.split(";")
        loops_clean = []
        for loop in loops_raw:
            if not loop:
                continue

            parts = loop.strip().split(":")
            if len(parts) < 3:
                continue

            start = parts[0]
            stop = parts[1]
            quality = parts[2]
            try:
                int(start)
                int(stop)
            except ValueError:
                continue

            loops_clean.append((quality, start, stop))

        if output_format == 'table':
            loops_clean = pd.DataFrame(loops_clean, columns=['Quality', 'Start_Position', 'End_Position'])

        return loops_clean
    
    def get_exons(self, ID, output_format = 'list'):
        """
        Return all exons associated with the ID in question, based on the primary APPRIS and/or MANE transcript associated with protein.


        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.

        Returns
        -------
        list of tuples, dict, int, or pandas.DataFrame
            tuples for each exon associated with the ID (exon_id, start_position_in_protein, end_position_in_protein, constitutive [1 = constitutive, 0 = alternative]). If output_format is 'table', returns a pandas DataFrame with columns ['Exon_Name', 'Start_Position', 'End_Position', 'Exon_ID'] instead of list.

            Returns an empty list if no exons are found


        Postconditions
        --------------
        Returns -1 if unable to find the ID

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        exons = record.get("exons", "")
        
        if not exons or exons.strip() == "":
            return []
                
        exons_raw = exons.split(";")
        exons_clean = []
        for i in exons_raw:
            if i:
                tmp = i.strip()
                parts = tmp.split(":")
                if len(parts) >= 3:
                    name = parts[0]
                    start = parts[1]
                    stop = parts[2]
                    constitutive = parts[3]
                    try:
                        float(start)
                        float(stop)
                    except ValueError:
                        continue
                    exons_clean.append((name, start, stop, constitutive))
                else:
                    print("ERROR: the exon structure did not match expected format %s"%(i))

        if output_format == 'table':
            exons_clean = pd.DataFrame(exons_clean, columns=['Exon_Name', 'Start_Position', 'End_Position', 'Constitutive'])
        return exons_clean
    
    def get_Scansite(self, ID):
        """
        DEPRECATED: Scansite predictions have been removed from the new data format.
        This method is kept for backward compatibility but will return an empty list.

        Returns
        -------
        empty list

        """
        return []

    def get_Scansite_byPos(self, ID, res_pos):
        """
        DEPRECATED: Scansite predictions have been removed from the new data format.
        This method is kept for backward compatibility but will return an empty dictionary.


        Returns
        -------
        empty dictionary

        """
        return {}

    def __parse_spyc_prediction_record(self, prediction_entry):
        """
        Parse one SpyC prediction entry in the format site:probability:class:confidence.

        Parameters
        ----------
        prediction_entry : str
            Single semicolon-delimited SpyC prediction entry.

        Returns
        -------
        tuple or None
            (site_position, probability, predicted_class, confidence) or None if
            the entry is malformed.
        """
        if not prediction_entry:
            return None

        parts = prediction_entry.strip().split(":")
        if len(parts) < 4:
            return None

        site_token, probability, predicted_class, confidence = parts[0], parts[1], parts[2], parts[3]

        site_match = re.search(r'\d+', str(site_token))
        if site_match is None:
            return None

        site_position = int(site_match.group())
        return (site_position, probability, predicted_class, confidence)

    def get_spyc_predictions(self, ID):
        """
        Return all SpyC predictions for a protein accession.

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        dict or int
            Dictionary keyed by site position with values
            [probability, class, confidence]. Returns -1 when SpyC predictions
            are absent for the protein or accession is not found.
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        raw_predictions = record.get("spyc_predictions", "")
        if raw_predictions is None or str(raw_predictions).strip() == "":
            return -1

        spyc_predictions = {}
        for prediction_entry in str(raw_predictions).split(";"):
            parsed_prediction = self.__parse_spyc_prediction_record(prediction_entry)
            if parsed_prediction is None:
                continue

            site_position, probability, predicted_class, confidence = parsed_prediction
            spyc_predictions[site_position] = [probability, predicted_class, confidence]

        if not spyc_predictions:
            return -1

        return spyc_predictions

    def get_spyc_predictions_byPos(self, ID, site):
        """
        Return SpyC prediction values for one site in a protein.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        site : int
            Residue position in the protein sequence.

        Returns
        -------
        tuple or int
            (probability, class, confidence) for the site. Returns -1 when no
            SpyC prediction exists for that site or accession.
        """
        try:
            site = int(site)
        except (TypeError, ValueError):
            return -1

        protein_predictions = self.get_spyc_predictions(ID)
        if protein_predictions == -1:
            return -1

        site_prediction = protein_predictions.get(site)
        if site_prediction is None:
            return -1

        return tuple(site_prediction)

    def _format_spyc_prediction_values(self, probability, confidence):
        """
        Convert raw SpyC probability/confidence into annotation-friendly values.

        Parameters
        ----------
        probability : object
            Raw probability value from SpyC predictions.
        confidence : object
            Raw confidence value from SpyC predictions.

        Returns
        -------
        tuple of str
            (confidence_value, probability_value) where probability_value uses
            "Training" for training-set entries (missing probability with
            present confidence).
        """
        confidence_value = "" if pd.isna(confidence) else str(confidence)
        if confidence_value.strip().lower() == 'nan':
            confidence_value = ""

        probability_is_missing = pd.isna(probability)
        if not probability_is_missing and str(probability).strip().lower() == 'nan':
            probability_is_missing = True

        if probability_is_missing:
            probability_value = "Training" if confidence_value else ""
        else:
            probability_value = str(probability)

        return confidence_value, probability_value

    def get_nearbyPTMs(self,ID,pos, window, output_format='list', include_hidden=False):
        """
        Return all PTMs within a specified window of a given position

        Parameters
        ----------
        ID : str
            SwissProt accession number
        pos : int
            Position in protein to search around
        window : int
            Number of residues upstream and downstream to include in the search
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        list of tuples or int
            tuple for each PTM within the specified window (position, residue, modification-type). If output_format is 'table', returns a pandas DataFrame with columns ['Position', 'Residue', 'Modification_Type'].
            Returns an empty list if no modifications are found within the window.

        Postconditions
        --------------
        If ID is not found, returns -1

        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        mods = self.get_PTMs(ID, include_hidden=include_hidden)
        nearby_mods = []
        for mod in mods:
            mod_pos = int(mod[0])
            if abs(mod_pos - pos) <= window:
                nearby_mods.append(mod)

        #convert to table if desired
        if output_format == 'table':
            # Convert to DataFrame
            nearby_df = pd.DataFrame(nearby_mods, columns=['Position', 'Residue', 'Modification_Type'])
            return nearby_df
        else:
            return nearby_mods

    def get_species(self,ID):
        """
        Return the species associated with the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        str
            Species associated with the ID

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return record["species"]
    
    def get_sequence(self, ID):
        """
        Return the sequence associated with the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        str
            Sequence associated with the ID
        
        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return record["sequence"]


    def get_gene_name(self, ID):
        """
        Return the gene_name of an ID
        
        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        tuple
            Gene name associated with the ID

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return (record["acc_gene"])
         

    
    def get_phosphosites(self,ID, output_format='list', include_hidden=False):
        """
        Return all phosphosites (S/T/Y phosphorylation) associated with the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        list of tuples or int
            tuple for each phosphosite associated with the ID (position, residue, modification-type). 
            Returns an empty list if no phosphosites are found.
        
        Postconditions
        --------------
        If ID is not found, returns -1
        """
        if output_format not in ['list', 'table']:
            raise ValueError("output_format must be 'list' or 'table'")
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        mods = self.get_PTMs(ID, output_format='list', include_hidden=include_hidden)
        phos_sites = []
        for mod in mods:
            residue = mod[1]
            mod_type = mod[2]
            if mod_type in ['Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine']:
                phos_sites.append(mod)
        #convert to table if desired
        if output_format == 'table':
            # Convert to DataFrame
            phos_df = pd.DataFrame(phos_sites, columns=['Position', 'Residue', 'Modification_Type'])
            return phos_df
        return phos_sites
    
    def get_region(self, position, regions):
        associated_region = regions[(regions['Start_Position'].astype(int) <= position) & (regions['End_Position'].astype(int) >= position)]
        return associated_region

    def _normalize_feature_records(self, feature_records):
        """
        Normalize feature tuples into (name, start, end, extra) form.
        """
        if feature_records in (-1, -2, None):
            return []

        normalized_records = []
        for feature in feature_records:
            if not feature or len(feature) < 3:
                continue

            name, start_position, end_position = feature[:3]
            extra = feature[3] if len(feature) > 3 else ""

            try:
                start_position = int(float(start_position))
                end_position = int(float(end_position))
            except (TypeError, ValueError):
                continue

            normalized_records.append((str(name), start_position, end_position, "" if extra is None else str(extra)))

        return normalized_records

    def _collect_site_annotation_context(self, ID):
        """
        Collect feature ranges required for site-level annotations.
        """
        return {
            'interpro': self._normalize_feature_records(self.get_domains(ID, domain_type='interpro', output_format='list')),
            'uniprot': self._normalize_feature_records(self.get_domains(ID, domain_type='uniprot', output_format='list')),
            'structure': self._normalize_feature_records(self.get_structure(ID, output_format='list')),
            'macro': self._normalize_feature_records(self.get_macro_molecular(ID, output_format='list')),
            'activation_loop': self._normalize_feature_records(self.get_activation_loops(ID, output_format='list')),
            'exons': self._normalize_feature_records(self.get_exons(ID, output_format='list')),
        }

    def _get_overlapping_features(self, features, position):
        """
        Return feature tuples that overlap a site position.
        """
        return [feature for feature in features if feature[1] <= position <= feature[2]]

    def _annotate_positions(self, ID, positions, context=None):
        """
        Build reusable site-level annotations for one accession and many positions.
        """
        if context is None:
            context = self._collect_site_annotation_context(ID)

        annotations = {
            'SpY-C Prediction': [],
            'SpY-C Prediction Probability': [],
        }
        for rule in self.POSITION_ANNOTATION_REGISTRY:
            annotations[rule['name_column']] = []
            if 'extra_column' in rule:
                annotations[rule['extra_column']] = []
            if rule.get('include_boolean'):
                annotations[rule['boolean_column']] = []

        def _append_empty_annotation():
            for rule in self.POSITION_ANNOTATION_REGISTRY:
                annotations[rule['name_column']].append("")
                if 'extra_column' in rule:
                    annotations[rule['extra_column']].append("")
                if rule.get('include_boolean'):
                    annotations[rule['boolean_column']].append(False)
            annotations['SpY-C Prediction'].append("")
            annotations['SpY-C Prediction Probability'].append("")

        for raw_position in positions:
            try:
                position = int(raw_position)
            except (TypeError, ValueError):
                _append_empty_annotation()
                continue

            for rule in self.POSITION_ANNOTATION_REGISTRY:
                hits = self._get_overlapping_features(context[rule['context_key']], position)
                annotations[rule['name_column']].append(';'.join(hit[0] for hit in hits))
                if 'extra_column' in rule:
                    annotations[rule['extra_column']].append(';'.join(hit[3] for hit in hits if hit[3]))
                if rule.get('include_boolean'):
                    annotations[rule['boolean_column']].append(bool(hits))

            spyc_prediction = self.get_spyc_predictions_byPos(ID, position)
            if spyc_prediction == -1:
                annotations['SpY-C Prediction'].append("")
                annotations['SpY-C Prediction Probability'].append("")
            else:
                probability, _, confidence = spyc_prediction
                confidence_value, probability_value = self._format_spyc_prediction_values(probability, confidence)
                annotations['SpY-C Prediction'].append(confidence_value)
                annotations['SpY-C Prediction Probability'].append(probability_value)

        return annotations

    def _join_non_empty(self, values):
        """
        Join list values, skipping empty/NaN-like entries.
        """
        cleaned_values = []
        for value in values:
            if value is None:
                continue
            value = str(value).strip()
            if not value or value.lower() == 'nan':
                continue
            cleaned_values.append(value)
        return ';'.join(cleaned_values)
    
    def get_annotated_PTMs(self, ID, include_hidden=False):
        """
        Given a UniProt ID, return a table of PTMs with annotations about whether they fall within domains, structures, or macro-molecular structures.
        
        Parameters
        ----------
        ID : str
            SwissProt accession number
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.
        
        Returns
        -------
        pd.DataFrame or int
            Dataframe with all PTMs associated with the ID, along with columns indicating domain names (InterPro and UniProt), structures, and macro-molecular structures that contain those PTMs

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        mods = self.get_PTMs(ID, output_format='table', include_hidden=include_hidden)
        if mods.empty:
            mods['evidence'] = pd.Series(dtype=object)
            mods['Domain_Names_InterPro'] = pd.Series(dtype=object)
            mods['InterPro_IDs'] = pd.Series(dtype=object)
            mods['Domain_Names_UniProt'] = pd.Series(dtype=object)
            mods['Structures'] = pd.Series(dtype=object)
            mods['Macro_Molecular_Structures'] = pd.Series(dtype=object)
            mods['In_Activation_Loop'] = pd.Series(dtype=bool)
            mods['Activation_Loop_Quality'] = pd.Series(dtype=object)
            mods['In_Exon'] = pd.Series(dtype=bool)
            mods['Exons'] = pd.Series(dtype=object)
            mods['Exon_Constitutive'] = pd.Series(dtype=object)
            mods['SpY-C Prediction'] = pd.Series(dtype=object)
            mods['SpY-C Prediction Probability'] = pd.Series(dtype=object)
            return mods

        _, filtered_evidence = self.__filter_modifications_by_visibility(record, include_hidden=include_hidden)
        mods['evidence'] = filtered_evidence

        site_annotations = self._annotate_positions(ID, mods['Position'].tolist())
        for column_name, values in site_annotations.items():
            mods[column_name] = values

        return mods



    def get_evidence(self, ID):
        """
        Return evidence scores associated with PTMs for the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        str
            Evidence scores associated with the ID
        
        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return record.get("evidence", "")

    # def get_mutations(self, ID):
    ## DEPRECATED: Mutations have been removed from the new data format.
    #     """
    #     Return all mutations associated with the ID
    #     """
    #     try:
    #         record = self.database[ID]
    #     except KeyError:
    #         return -1
        
    #     mutations = record.get("mutations", "")
    #     if not mutations:
    #         return []
        
    #     mut_raw = mutations.split(";")
    #     mut_clean = []
    #     for i in mut_raw:
    #         if i.strip():
    #             mut_clean.append(i.strip())
    #     return mut_clean

    def get_GO(self, ID):
        """
        Return all GO terms associated with the ID in question

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        list of tuples or int
            tuple for each GO term associated with the ID (GO_term, type). type is 'F' (Molecular Function), 'P' (Biological Process), or 'C' (Cellular Component). 
            Returns an empty list if no GO terms are found.

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        GO_terms = record["GO_terms"]
        if len(GO_terms)==0:
            return []

        GO_termsArr = GO_terms.split(";")
        GO_terms_clean = []
        for i in GO_termsArr:
            term = i.strip()
            if term:
                # Extract GO type (last character should be F, P, or C)
                if len(term) > 0 and term[-1] in ['F', 'P', 'C']:
                    go_type = term[-1]
                    go_id = term[:-1].strip()
                    GO_terms_clean.append((go_id, go_type))
                else:
                    # Fallback if format doesn't match expected
                    GO_terms_clean.append((term, ''))

        return GO_terms_clean
    
    def get_all_protein_info(self, ID):
        """
        Return all available information for the ID as a dictionary

        Parameters
        ----------
        ID : str
            SwissProt accession number

        Returns
        -------
        dict or int
            Dictionary containing all available information for the ID.

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        return record

    # def get_kinaseLoops(self, ID):
    # DEPRECATED: Kinase loop information has been removed from the new data format.
    #     """
    #     Return kinase loop information for the ID
    #     """
    #     try:
    #         record = self.database[ID]
    #     except KeyError:
    #         return -1
    #     return record.get("kinase_loops", "")

    def get_accessions(self,ID):
        """
        Return all accession numbers associated with the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number
        
        Returns
        -------
        list of str or int
            List of accession IDs associated with the Swissprot accession.

        Postconditions
        --------------
        If unable to find the ID, returns -1


        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        accessions = record["accessions"].split(";")
        return [acc.strip() for acc in accessions]

    # def get_PTMs_withEvidenceThreshold(self, ID, evidenceThreshold):
    #     """
    #     Return PTMs with evidence scores above a threshold
    #     """
    #     try:
    #         record = self.database[ID]
    #     except KeyError:
    #         return -1
        
    #     # Implementation depends on your evidence format
    #     return self.get_PTMs(ID)

    def get_PTMs_withEvidence(self, ID, include_hidden=False):
        """
        Return PTMs with their associated evidence information

        Parameters
        ----------
        ID : str
            SwissProt accession number
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        tuple or int
            Tuple containing (mods, evidence) where mods is a list of PTMs and evidence is the associated evidence information.

        Postconditions
        --------------
        If unable to find the ID, returns -1
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        mods = self.get_PTMs(ID, include_hidden=include_hidden)
        _, filtered_evidence = self.__filter_modifications_by_visibility(record, include_hidden=include_hidden)
        evidence = ';'.join(filtered_evidence)
        return (mods, evidence)
    
    def return_species_nr_uniprot_ids(self):
        """
        Return a dictionary with species as keys and number of unique uniprot IDs as values. 
        Keep only those that are flagged as uniprot IDs in in the non-redundant list.

        Returns
        -------
        species_dict : dict
            Dictionary with species names as keys and lists of uniprot IDs as values. {species_name: [list of uniprot IDs]}
        species_reference_bool : dict
            Dictionary indicating whether each species is a reference species and contains all protein IDs (True/False). If False, it means that only nonredundant uniprot records with PTMs are included in the list.

        """
        species_dict = {}
        for ID in self.database.keys():
            species = self.get_species(ID)
            non_redundant = self.database[ID].get("swissprot_nr").strip()
            if non_redundant:
                if species not in species_dict:
                    species_dict[species] = set()
                species_dict[species].add(ID)
        
        # at the end, convert sets to lists 
        for species in species_dict:
            species_dict[species] = list(species_dict[species])

        species_reference_bool = {
                        "Homo sapiens": True,
                        "Mus musculus": True,
                        "Rattus norvegicus": False,
                        "Bos taurus": False,
                        "Drosophila melanogaster"  : False,
                        "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)": False,
                    }
        
        return species_dict, species_reference_bool

    def search_by_peptide(self, peptide, include_hidden=False):
        """
        Search the ProteomeScout dataset for proteins containing a specific peptide sequence.

        Parameters
        ----------
        peptide : str
            Peptide sequence to search for (case-insensitive). Modification indicators (lowercased residues) are ignored for matching.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments when
            counting matched protein PTMs. Defaults to False.

        Returns
        -------
        tuple
            A tuple containing:
            - list of str: UniProt accessions (IDs) of proteins containing the peptide
            - dict: Dictionary with accessions as keys and protein information as values. Each value contains:
                - 'species': Species of the protein
                - 'num_ptms': Number of PTMs (modifications) in the protein
                - 'gene_name': Gene name associated with the protein
                - 'protein_name': Name of the protein
                - 'protein_sequence': Full protein sequence

        Notes
        -----
        The search is case-insensitive and removes modification notation (lowercased residues) from the input peptide before searching.
        For example, a peptide "sPEPTIDE" will be converted to "SPEPTIDE" before searching.

        Examples
        --------
        >>> accessions, info_dict = api.search_by_peptide('PEPTIDE')
        >>> print(accessions)
        ['P12345', 'P67890']
        >>> print(info_dict['P12345'])
        {'species': 'Homo sapiens', 'num_ptms': 5, 'gene_name': 'GENE1', 'protein_name': 'Protein One', 'protein_sequence': 'MXXXXX...'}
        """
        # Normalize the search peptide: convert to uppercase for case-insensitive search
        search_peptide = peptide.upper()
        
        # Remove modification indicators (lowercased residues in the original peptide)
        # by keeping only uppercase letters
        search_peptide_clean = ''.join([c for c in search_peptide if c.isalpha()])
        
        matching_accessions = []
        matching_info = {}
        
        # Search through all proteins in the database
        for accession in self.uniqueKeys:
            record = self.database[accession]
            
            # Get the protein sequence and convert to uppercase for comparison
            sequence = record.get("sequence", "").upper()
            
            # Check if the peptide is in the sequence
            if search_peptide_clean in sequence:
                matching_accessions.append(accession)
                
                # Get PTM count
                mods_list = self.get_PTMs(accession, output_format='list', include_hidden=include_hidden)
                if mods_list not in (-1, []) and mods_list:
                    num_ptms = len(mods_list)
                else:
                    num_ptms = 0
                
                # Compile information for this accession
                matching_info[accession] = {
                    'species': record.get("species", ""),
                    'num_ptms': num_ptms,
                    'gene_name': record.get("acc_gene", ""),
                    'protein_name': record.get("protein_name", ""),
                    'protein_sequence': record.get("sequence", "")
                }
        
        return matching_accessions, matching_info


class SpeciesReferenceDataset(ProteomeScoutAPI):
    """
    Build PTM-centric reference datasets for nonredundant ProteomeScout species.

    Each row in the output represents one PTM site on one UniProt protein and
    includes sequence context plus feature overlap annotations.
    """

    OUTPUT_COLUMNS = [
        'species',
        'gene_name',
        'uniprot_id',
        'site',
        'ptm_type',
        'oriented_peptide',
        'in_interpro_domain',
        'interpro_domains',
        'in_structure',
        'structures',
        'in_macro_molecular_structure',
        'macro_molecular_structures',
        'in_activation_loop',
        'activation_loop_qualities',
        'in_exon',
        'exons',
        'exon_constitutive',
        'SpY-C Prediction',
        'SpY-C Prediction Probability',
    ]

    OPTIONAL_ANNOTATION_COLUMNS = [
        'SpY-C Prediction',
        'SpY-C Prediction Probability',
        'in_exon',
        'exons',
        'exon_constitutive',
    ]

    def __init__(self, flank=7, version=config.VERSION, update=config.UPDATE):
        super().__init__(version=version, update=update)
        self.flank = flank

    def _build_oriented_peptide(self, sequence, position, residue):
        if not sequence:
            return ""

        sequence = sequence.upper()
        if position < 1 or position > len(sequence):
            return ""

        zero_based_position = position - 1
        left = sequence[max(0, zero_based_position - self.flank):zero_based_position]
        right = sequence[zero_based_position + 1:zero_based_position + self.flank + 1]

        left = "_" * (self.flank - len(left)) + left
        right = right + "_" * (self.flank - len(right))

        center_residue = sequence[zero_based_position]
        if residue and center_residue != str(residue).upper():
            center_residue = str(residue).upper()

        return f"{left}{center_residue.lower()}{right}"

    def _species_to_filename(self, species):
        safe_species = "".join(
            character.lower() if character.isalnum() else "_"
            for character in species
        )
        while "__" in safe_species:
            safe_species = safe_species.replace("__", "_")
        return safe_species.strip("_") + "_reference_dataset.csv"

    def _drop_empty_optional_annotation_columns(self, reference_table):
        """
        Remove optional annotation columns that contain no annotations.

        Parameters
        ----------
        reference_table : pd.DataFrame
            Species-level PTM-centric reference table.

        Returns
        -------
        pd.DataFrame
            Table with empty optional annotation columns removed.
        """
        if reference_table is None or reference_table.empty:
            return reference_table

        columns_to_drop = []
        for column_name in self.OPTIONAL_ANNOTATION_COLUMNS:
            if column_name not in reference_table.columns:
                continue

            column_values = reference_table[column_name]
            if pd.api.types.is_bool_dtype(column_values):
                has_annotations = column_values.fillna(False).any()
            else:
                cleaned_values = column_values.fillna('').astype(str).str.strip()
                has_annotations = (cleaned_values != '').any() and (~cleaned_values.str.lower().eq('nan')).any()

            if not has_annotations:
                columns_to_drop.append(column_name)

        if columns_to_drop:
            reference_table = reference_table.drop(columns=columns_to_drop)

        return reference_table

    def build_protein_reference_dataset(self, uniprot_id, include_hidden=False):
        """
        Build a PTM-centric reference table for a single UniProt accession.

        Parameters
        ----------
        uniprot_id : str
            UniProt accession to convert into PTM-centric rows.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        pd.DataFrame or int
            DataFrame with one row per PTM site. Returns -1 if the accession is
            not found.
        """
        try:
            record = self.database[uniprot_id]
        except KeyError:
            return -1

        annotated_ptms = self.get_annotated_PTMs(uniprot_id, include_hidden=include_hidden)
        if isinstance(annotated_ptms, int) and annotated_ptms == -1:
            return -1
        if annotated_ptms.empty:
            return pd.DataFrame(columns=self.OUTPUT_COLUMNS)

        sequence = self.get_sequence(uniprot_id)
        gene_name = self.get_gene_name(uniprot_id)
        species = record.get('species', '')

        def _value_or_empty(row, key):
            value = row.get(key, "")
            if pd.isna(value):
                return ""
            return value

        reference_rows = []
        for _, ptm_row in annotated_ptms.iterrows():
            try:
                site_position = int(ptm_row['Position'])
            except (TypeError, ValueError):
                continue

            residue = _value_or_empty(ptm_row, 'Residue')
            ptm_type = _value_or_empty(ptm_row, 'Modification_Type')
            interpro_hits = _value_or_empty(ptm_row, 'Domain_Names_InterPro')
            structure_hits = _value_or_empty(ptm_row, 'Structures')
            macro_hits = _value_or_empty(ptm_row, 'Macro_Molecular_Structures')
            activation_loop_hits = _value_or_empty(ptm_row, 'Activation_Loop_Quality')
            exon_hits = _value_or_empty(ptm_row, 'Exons')
            exon_constitutive = _value_or_empty(ptm_row, 'Exon_Constitutive')
            spyc_confidence = _value_or_empty(ptm_row, 'SpY-C Prediction')
            spyc_probability = _value_or_empty(ptm_row, 'SpY-C Prediction Probability')

            in_activation_loop = bool(ptm_row.get('In_Activation_Loop', False))
            in_exon = bool(ptm_row.get('In_Exon', False))

            reference_rows.append({
                'species': species,
                'gene_name': gene_name,
                'uniprot_id': uniprot_id,
                'site': f"{residue}{site_position}",
                'ptm_type': ptm_type,
                'oriented_peptide': self._build_oriented_peptide(sequence, site_position, residue),
                'in_interpro_domain': bool(interpro_hits),
                'interpro_domains': interpro_hits,
                'in_structure': bool(structure_hits),
                'structures': structure_hits,
                'in_macro_molecular_structure': bool(macro_hits),
                'macro_molecular_structures': macro_hits,
                'in_activation_loop': in_activation_loop,
                'activation_loop_qualities': activation_loop_hits,
                'in_exon': in_exon,
                'exons': exon_hits,
                'exon_constitutive': exon_constitutive,
                'SpY-C Prediction': spyc_confidence,
                'SpY-C Prediction Probability': spyc_probability,
            })

        return pd.DataFrame(reference_rows, columns=self.OUTPUT_COLUMNS)

    def build_species_reference_dataset(self, species, output_file=None, include_hidden=False):
        """
        Build a PTM-centric reference table for a nonredundant species.

        Parameters
        ----------
        species : str
            Species name present in return_species_nr_uniprot_ids().
        output_file : str, optional
            If provided, write the resulting DataFrame to CSV.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        pd.DataFrame
            PTM-centric reference table for the requested species.
        """
        species_dict, _ = self.return_species_nr_uniprot_ids()
        if species not in species_dict:
            raise KeyError(f"Species '{species}' not found in nonredundant species list")

        protein_reference_tables = []
        for uniprot_id in species_dict[species]:
            protein_reference_table = self.build_protein_reference_dataset(uniprot_id, include_hidden=include_hidden)
            if isinstance(protein_reference_table, pd.DataFrame) and not protein_reference_table.empty:
                protein_reference_tables.append(protein_reference_table)

        if protein_reference_tables:
            species_reference_table = pd.concat(protein_reference_tables, ignore_index=True)
        else:
            species_reference_table = pd.DataFrame(columns=self.OUTPUT_COLUMNS)

        species_reference_table = self._drop_empty_optional_annotation_columns(species_reference_table)

        if output_file is not None:
            species_reference_table.to_csv(output_file, index=False)

        return species_reference_table

    def write_all_species_reference_datasets(self, output_dir, include_hidden=False):
        """
        Write one PTM-centric reference CSV per nonredundant species.

        Parameters
        ----------
        output_dir : str
            Directory where species-specific CSV files should be written.
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments.
            Defaults to False.

        Returns
        -------
        dict
            Mapping of species name to the written CSV path.
        """
        os.makedirs(output_dir, exist_ok=True)

        species_dict, _ = self.return_species_nr_uniprot_ids()
        output_files = {}
        for species in sorted(species_dict):
            output_path = os.path.join(output_dir, self._species_to_filename(species))
            self.build_species_reference_dataset(species, output_file=output_path, include_hidden=include_hidden)
            output_files[species] = output_path

        return output_files

    
class ProteomicDataset(ProteomeScoutAPI):
    """
    Class for annotating phosphoproteomic datasets with gene-level and site-specific information. Inherits methods from ProteomeScoutAPI.

    Parameters
    ----------
    dataset: pd.DataFrame
        phosphoproteomic dataset to annotate
    accession_col: str
        column name containing SwissProt accessions
    peptide_col: str
        column name containing formatted peptides, with modification sites lowercased
    find_site: bool
        whether to find modification sites within protein (based on which residues are lowercased in peptide)
    domain_source: str
        source of domain annotations ('interpro' or 'uniprot')
    GO_terms: bool
        whether to include GO term annotations
    version : int, optional
        Version number of the ProteomeScout dataset to use.
    update : bool, optional
        Whether to update the local dataset when a version mismatch is detected.
    
    
    """
    def __init__(self, dataset, accession_col = 'acc', peptide_col = 'pep', find_site = True, domain_source = 'interpro', GO_terms = True, version=config.VERSION, update=config.UPDATE):
        #check to make sure columns exist
        if accession_col not in dataset.columns:
            raise KeyError(f"Accession column '{accession_col}' not found in DataFrame")
        if peptide_col not in dataset.columns:
            raise KeyError(f"Peptide column '{peptide_col}' not found in DataFrame")
        
        super().__init__(version=version, update=update)
        self.dataset = dataset
        self.accession_col = accession_col
        self.peptide_col = peptide_col
        self.find_site = find_site
        self.domain_source = domain_source
        self.GO_terms = GO_terms

    def check_phosphosites(self, accessions, positions, include_hidden=False):
        """
        Check if positions are documented phosphosites in ProteomeScout for the given accessions
        
        Parameters
        ----------
        accessions : str
            SwissProt accession number
        positions : list of int
            List of positions to check
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments when
            checking documented phosphosites. Defaults to False.

        Returns
        -------
        str
            Semicolon-separated string of 1s and 0s indicating whether each position is a documented phosphosite (1) or not (0)
        """
        #now check if it's been annotated as a known modification site before
        phosphosites = self.get_phosphosites(accessions, include_hidden=include_hidden)
        
        found_arr = []
        for pos in positions:
            found = 0
            for mod in phosphosites:
                (pos_mod, res_mod, type_mod) = mod
                #print("Debug, comparing %d to %d"%(pos, int(pos_mod)))
                if pos == int(pos_mod):
                    found = 1
                    #print("FOUND IT!")
                    break

            found_arr.append(str(found))
        #combine into string
        documented_sites = ';'.join(found_arr)
        return documented_sites
    
    def get_domains_with_site(self, domains, positions, domain_descriptor = 'name'):
        """
        Check if positions are within any domains

        Parameters
        ----------
        domains : list of tuples
            List of domains (domain_name, start_position, end_position, domain_id)
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of domain names that contain the positions
        """
        in_domain_arr = []
        for pos in positions:
            for domain in domains:
                domain_name, domain_start, domain_stop, domain_id = domain
                if pos >= int(domain_start) and pos <= int(domain_stop):
                    if domain_descriptor == 'name':
                        in_domain_arr.append(domain_name)
                    elif domain_descriptor == 'id':
                        in_domain_arr.append(domain_id)
                    else:
                        raise ValueError("domain_descriptor must be 'name' or 'id'")
                    continue #continue back to the next position

        site_in_domain = ';'.join(in_domain_arr)
        return site_in_domain
    
    def get_macro_with_site(self, macro_mol, positions):
        """
        Check if positions are within any macro-molecular structures

        Parameters
        ----------
        macro_mol : list of tuples
            List of macro-molecular structures (macro_name, start_position, end_position)
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of macro-molecular structure names that contain the positions
        """
        in_macro_arr = []
        for pos in positions:
            for macro in macro_mol:
                macro_name, macro_start, macro_stop = macro
                if pos >= int(macro_start) and pos <= int(macro_stop):
                    in_macro_arr.append(macro_name)
                    continue #continue back to the next position

        site_in_macro = ';'.join(in_macro_arr)
        return site_in_macro
    
    def get_structure_with_site(self, structure, positions):
        """
        Check if positions are within any macro-molecular structures

        Parameters
        ----------
        macro_mol : list of tuples
            List of macro-molecular structures (macro_name, start_position, end_position)
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of macro-molecular structure names that contain the positions
        """
        in_macro_arr = []
        for pos in positions:
            for macro in structure:
                macro_name, macro_start, macro_stop = macro
                if pos >= int(macro_start) and pos <= int(macro_stop):
                    in_macro_arr.append(macro_name)
                    continue #continue back to the next position

        site_in_macro = ';'.join(in_macro_arr)
        return site_in_macro

    def get_activation_loops_with_site(self, activation_loops, positions):
        """
        Check if positions are within any activation loop regions.

        Parameters
        ----------
        activation_loops : list of tuples
            List of activation loops (quality, start_position, end_position)
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of activation loop qualities for matching
            positions.
        """
        in_activation_loop_arr = []
        for pos in positions:
            for loop in activation_loops:
                quality, loop_start, loop_stop = loop
                if pos >= int(loop_start) and pos <= int(loop_stop):
                    in_activation_loop_arr.append(quality)
                    continue

        site_in_activation_loop = ';'.join(in_activation_loop_arr)
        return site_in_activation_loop
    
    def get_exons_with_site(self, exons, positions):
        """
        Check if positions are within any exon regions.

        Parameters
        ----------
        exons : list of tuples
            List of exons (exon_id, start_position, end_position, constitutive)
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of exon IDs for matching positions.
        """
        in_exon_arr = []
        for pos in positions:
            for exon in exons:
                exon_id, exon_start, exon_stop, constitutive = exon
                if pos >= float(exon_start) and pos <= float(exon_stop):
                    in_exon_arr.append(exon_id)
                    continue

        site_in_exon = ';'.join(in_exon_arr)
        return site_in_exon


    def annotate_peptide(self, accession, peptide, include_hidden=False):
        """
        Given a SwissProt accession and peptide sequence, annotate with gene-level and site-specific information from ProteomeScout.

        Parameters
        ----------
        accession : str
            SwissProt accession number
        peptide : str
            Peptide sequence with modification sites lowercased
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments when
            checking documented phosphosites. Defaults to False.

        Returns
        -------
        dict or int
            Dictionary containing annotation information:
            - 'gene_name': Gene name associated with the protein accession.
            - 'domains': Semicolon-separated string of domain names associated with the protein.
            - 'domain_architecture': String representation of the domain architecture (order of domains)
            - 'exons': Semicolon-separated string of exon IDs associated with the protein.
            - 'GO_terms': Semicolon-separated string of GO terms associated with the protein.

            If find_site is True, additional keys are included:
            - 'modification_sites': Semicolon-separated string of modification sites found in the peptide.
            - 'aligned_peps': Aligned peptide sequences found in the protein sequence (if find_site is True).
            - 'documented_phosphosites': Semicolon-separated string indicating whether each modification site is documented (1) or not (0).
            - 'site_in_domain': Semicolon-separated string of domain names that contain the modification sites
            - 'site_in_macro': Semicolon-separated string of macro-molecular structure names that contain the modification sites
            - 'site_in_structure': Semicolon-separated string of structure names that contain the modification sites
            - 'site_in_activation_loop': Semicolon-separated string of activation loop qualities for matching positions
            - 'site_in_exon': Semicolon-separated string of exon IDs for matching positions
            - 'SpY-C Prediction': Semicolon-separated SpyC confidence values at matching modification sites.
            - 'SpY-C Prediction Probability': Semicolon-separated SpyC probability values at matching sites. Returns "Training" when a site is from the training set.
            Returns -1 if unable to find the accession in the database.
        """
        seq = self.get_sequence(accession)
        gene_name = self.get_gene_name(accession)

        if seq == -1: 
            #will add -1 to information returned
            return -1

        domains = self.get_domains(accession, domain_type = self.domain_source)
        #DOMAINS
        #get domains and the domain strings
        domain_string, domain_architecture = helpers.returnDomainArchString(domains)

        #GO TERMS, get them as a string
        GO_terms = self.database[accession]['GO_terms']

        if self.find_site:
            #based on peptide, find position(s) in sequence (also suppress stdout from helper function)
            with contextlib.redirect_stdout(None):
                alignedPeps, seqPosArr, aaArr = helpers.returnOrientedPhosphoPeptide(seq, peptide)

            if len(alignedPeps) == 0:
                mod_sites = np.nan
                aligned_peptides = np.nan
                documented_sites = np.nan
                sites_in_domain = np.nan
                sites_in_domain_id = np.nan
                sites_in_macro = np.nan
                sites_in_structure = np.nan
                sites_in_activation_loop = np.nan
                sites_in_exon = np.nan
                site_spyc_prediction = np.nan
                site_spyc_probability = np.nan
            else:
                pos_aa_arr = []
                for i in range(0, len(seqPosArr)):
                    pos_aa_arr.append(aaArr[i]+str(seqPosArr[i]))

                mod_sites = ';'.join(pos_aa_arr)
                aligned_peptides = ';'.join(alignedPeps)

                site_annotations = self._annotate_positions(accession, seqPosArr)

                #now for each modification, ask if it's in a domain:
                sites_in_domain = self._join_non_empty(site_annotations['Domain_Names_InterPro'])
                sites_in_domain_id = self._join_non_empty(site_annotations['InterPro_IDs'])

                #or if it's in a macro-molecular structure
                sites_in_macro = self._join_non_empty(site_annotations['Macro_Molecular_Structures'])

                #or if it's in a structure
                sites_in_structure = self._join_non_empty(site_annotations['Structures'])

                #or if it's in an activation loop
                sites_in_activation_loop = self._join_non_empty(site_annotations['Activation_Loop_Quality'])

                #and the primary exon it is located in
                sites_in_exon = self._join_non_empty(site_annotations['Exons'])

                site_spyc_prediction = ';'.join(site_annotations['SpY-C Prediction'])
                site_spyc_probability = ';'.join(site_annotations['SpY-C Prediction Probability'])


                #now check if it's been annotated as a known modification site before
                documented_sites = self.check_phosphosites(accession, seqPosArr, include_hidden=include_hidden)



            output = {
                'gene_name': gene_name,
                'domains': domain_string,
                'domain_architecture': domain_architecture,
                'GO_terms': GO_terms,
                'aligned_peps': aligned_peptides,
                'modification_sites': mod_sites,
                'documented_phosphosites': documented_sites,
                'site_in_domain:name': sites_in_domain,
                'site_in_domain:interpro': sites_in_domain_id,
                'site_in_macro': sites_in_macro,
                'site_in_structure': sites_in_structure,
                'site_in_activation_loop': sites_in_activation_loop,
                'site_in_exon': sites_in_exon,
                'SpY-C Prediction': site_spyc_prediction,
                'SpY-C Prediction Probability': site_spyc_probability,
                }
        else:
            output = {
                'gene_name': gene_name,
                'domains': domain_string,
                'domain_architecture': domain_architecture,
                'GO_terms': GO_terms
                }
            
        return output
    
    def annotate_dataset(self, include_hidden=False):
        """
        Given a proteomic dataset in self.dataset, annotate each row with information from ProteomeScout.
        
        Postconditions
        --------------
        The dataset DataFrame is updated in place with new columns containing annotation information. This includes:
        - 'gene_name': Gene name associated with the protein accession.
        - 'domains': Semicolon-separated string of domain names associated with the protein.
        - 'domain_architecture': String representation of the domain architecture (order of domains)
        - 'GO_terms': Semicolon-separated string of GO terms associated with the protein.

        If find_site is True, additional columns are added:
        - 'modification_sites': Semicolon-separated string of modification sites found in the peptide.
        - 'aligned_peps': Aligned peptide sequences found in the protein sequence (if find_site is True).
        - 'documented_phosphosites': Semicolon-separated string indicating whether each modification site is documented (1) or not (0).
        - 'site_in_domain': Semicolon-separated string of domain names that contain the modification sites
        - 'site_in_macro': Semicolon-separated string of macro-molecular structure names that contain the modification sites
        - 'site_in_structure': Semicolon-separated string of structure names that contain the modification sites
        - 'site_in_activation_loop': Semicolon-separated string of activation loop qualities for matching positions
        include_hidden : bool, optional
            Whether to include PTMs supported only by hidden experiments when
            checking documented phosphosites. Defaults to False.
        """
        new_info = []
        for i, row in self.dataset.iterrows():
            acc = row[self.accession_col]
            pep = row[self.peptide_col]
            annotation = self.annotate_peptide(acc, pep, include_hidden=include_hidden)
            if annotation == -1:
                #if sequence not found, skip
                new_info.append(pd.Series())
            else:
                #if found, convert to series and add to list
                annotation = pd.Series(annotation, name=i)
                if self.find_site:
                    if annotation['aligned_peps'] != annotation['aligned_peps']:  # Check for NaN
                        annotation['PScout_Errors'] = f"Peptide not found in sequence"
                    else:
                        annotation['PScout_Errors'] = np.nan
                new_info.append(annotation)

        #combine row information and add to dataframe
        new_info_df = pd.DataFrame(new_info)
        self.dataset = pd.concat([self.dataset, new_info_df], axis=1)
        self.dataset['PScout_Errors'] = self.dataset.apply(lambda row: row['PScout_Errors'] if row['gene_name'] == row['gene_name'] else 'Accession not found in database', axis = 1)

        print("Annotation information has been appended to the dataset DataFrame as new columns.")