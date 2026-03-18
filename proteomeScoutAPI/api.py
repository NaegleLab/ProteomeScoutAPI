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
import os, json, contextlib
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
                first_line = f.readline()
                
            if not len(first_line.split("\t")) == 17:
                raise BadProteomeScoutFile("N/A")
            
                
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
        headers_raw = content[0].split('\t')

        headers=[]
        for i in headers_raw:
            headers.append(i.strip())

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
            uniprot_id = record[12].strip() if len(record) > 12 and record[12].strip() else ProteomeScoutID
                
            # add to the list of unique keys
            self.uniqueKeys.append(uniprot_id)

            # now construct the object dictionary
            OBJ={}
            for i in range(1,17):
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

    def get_PTMs(self, ID, output_format='list'):
        """
        Return all PTMs associated with the ID in question.

        Parameters
        ----------
        ID : str
            SwissProt accession number
        output_format : str, optional
            Format of the output ('list' or 'table'). Defaults to 'list'.

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

        mods = record["modifications"]
        if not mods or mods.strip() == "":
            return []
        
        mods_clean = helpers.clean_PTM_string(mods)

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

    def get_nearbyPTMs(self,ID,pos, window, output_format='list'):
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
        
        mods = self.get_PTMs(ID)
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
         

    
    def get_phosphosites(self,ID, output_format='list'):
        """
        Return all phosphosites (S/T/Y phosphorylation) associated with the ID

        Parameters
        ----------
        ID : str
            SwissProt accession number

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
        
        mods = self.get_PTMs(ID, output_format='list')
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
    
    def get_annotated_PTMs(self, ID):
        """
        Given a UniProt ID, return a table of PTMs with annotations about whether they fall within domains, structures, or macro-molecular structures.
        
        Parameters
        ----------
        ID : str
            SwissProt accession number
        
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
        
        mods = self.get_PTMs(ID, output_format='table')
        mods['evidence'] = self.get_evidence(ID).split(';')

        #extract domains
        domains = self.get_domains(ID, output_format='table')
        domains['interpro'] = domains['interpro'].astype({'Start_Position': int, 'End_Position': int})
        domains['uniprot'] = domains['uniprot'].astype({'Start_Position': int, 'End_Position': int})

        #extract structures
        structure = self.get_structure(ID, output_format='table')
        structure = structure.astype({'Start_Position': int, 'End_Position': int})

        #extract macro-molecular structures
        macro_mol = self.get_macro_molecular(ID, output_format='table')
        macro_mol = macro_mol.astype({'Start_Position': int, 'End_Position': int})

        site_domain_name = {'interpro': [], 'uniprot': []}
        site_interpro_ids = []
        site_structures = []
        site_macro = []
        for i, row in mods.iterrows():
            pos = int(row['Position'])
            # Check if in domain
            for dbase in domains.keys():
                trim_domains = domains[dbase]
                trim_domains = trim_domains.loc[(trim_domains['Start_Position'] <= pos) & (trim_domains['End_Position'] >= pos)].copy()
                if not trim_domains.empty:
                    site_domain_name[dbase].append(";".join(trim_domains['Domain_Name'].tolist()))
                    if dbase == 'interpro':
                        site_interpro_ids.append(";".join(trim_domains['Domain_ID'].tolist()))
                else:
                    site_domain_name[dbase].append("")
                    if dbase == 'interpro':
                        site_interpro_ids.append("")



            # Check if in structure
            trim_structure = structure.loc[(structure['Start_Position'] <= pos) & (structure['End_Position'] >= pos)].copy()
            if not trim_structure.empty:
                site_structures.append(";".join(trim_structure['Structure_Name'].tolist()))
            else:
                site_structures.append("")

            # Check if in macro-molecular structure
            trim_macro = macro_mol.loc[(macro_mol['Start_Position'] <= pos) & (macro_mol['End_Position'] >= pos)].copy()
            if not trim_macro.empty:
                site_macro.append(";".join(trim_macro['Macro_Name'].tolist()))
            else:
                site_macro.append("")

        mods['Domain_Names_InterPro'] = site_domain_name['interpro']
        mods['InterPro_IDs'] = site_interpro_ids
        mods['Domain_Names_UniProt'] = site_domain_name['uniprot']
        mods['Structures'] = site_structures
        mods['Macro_Molecular_Structures'] = site_macro
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

    def get_PTMs_withEvidence(self, ID):
        """
        Return PTMs with their associated evidence information

        Parameters
        ----------
        ID : str
            SwissProt accession number

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
        
        mods = self.get_PTMs(ID)
        evidence = self.get_evidence(ID)
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
    
    
    """
    def __init__(self, dataset, accession_col = 'acc', peptide_col = 'pep', find_site = True, domain_source = 'interpro', GO_terms = True):
        #check to make sure columns exist
        if accession_col not in dataset.columns:
            raise KeyError(f"Accession column '{accession_col}' not found in DataFrame")
        if peptide_col not in dataset.columns:
            raise KeyError(f"Peptide column '{peptide_col}' not found in DataFrame")
        
        super().__init__()
        self.dataset = dataset
        self.accession_col = accession_col
        self.peptide_col = peptide_col
        self.find_site = find_site
        self.domain_source = domain_source
        self.GO_terms = GO_terms

    def check_phosphosites(self, accessions, positions):
        """
        Check if positions are documented phosphosites in ProteomeScout for the given accessions
        
        Parameters
        ----------
        accessions : str
            SwissProt accession number
        positions : list of int
            List of positions to check

        Returns
        -------
        str
            Semicolon-separated string of 1s and 0s indicating whether each position is a documented phosphosite (1) or not (0)
        """
        #now check if it's been annotated as a known modification site before
        phosphosites = self.get_phosphosites(accessions)
        
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

    def annotate_peptide(self, accession, peptide):
        """
        Given a SwissProt accession and peptide sequence, annotate with gene-level and site-specific information from ProteomeScout.

        Parameters
        ----------
        accession : str
            SwissProt accession number
        peptide : str
            Peptide sequence with modification sites lowercased

        Returns
        -------
        dict or int
            Dictionary containing annotation information:
            - 'gene_name': Gene name associated with the protein accession.
            - 'domains': Semicolon-separated string of domain names associated with the protein.
            - 'domain_architecture': String representation of the domain architecture (order of domains)
            - 'GO_terms': Semicolon-separated string of GO terms associated with the protein.

            If find_site is True, additional keys are included:
            - 'modification_sites': Semicolon-separated string of modification sites found in the peptide.
            - 'aligned_peps': Aligned peptide sequences found in the protein sequence (if find_site is True).
            - 'documented_phosphosites': Semicolon-separated string indicating whether each modification site is documented (1) or not (0).
            - 'site_in_domain': Semicolon-separated string of domain names that contain the modification sites
            - 'site_in_macro': Semicolon-separated string of macro-molecular structure names that contain the modification sites

            Returns -1 if unable to find the accession in the database.
        """
        seq = self.get_sequence(accession)
        gene_name = self.get_gene_name(accession)

        if seq == -1: 
            #will add -1 to information returned
            return -1

        domains = self.get_domains(accession, domain_type = self.domain_source)
        macro = self.get_macro_molecular(accession)

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
            else:
                pos_aa_arr = []
                for i in range(0, len(seqPosArr)):
                    pos_aa_arr.append(aaArr[i]+str(seqPosArr[i]))

                mod_sites = ';'.join(pos_aa_arr)
                aligned_peptides = ';'.join(alignedPeps)

                #now for each modification, ask if it's in a domain:
                sites_in_domain = self.get_domains_with_site(domains, seqPosArr)
                sites_in_domain_id = self.get_domains_with_site(domains, seqPosArr, domain_descriptor='id')

                #or if it's in a macro-molecular structure
                sites_in_macro = self.get_macro_with_site(macro, seqPosArr)


                #now check if it's been annotated as a known modification site before
                documented_sites = self.check_phosphosites(accession, seqPosArr)



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
                'site_in_macro': sites_in_macro
                }
        else:
            output = {
                'gene_name': gene_name,
                'domains': domain_string,
                'domain_architecture': domain_architecture,
                'GO_terms': GO_terms
                }
            
        return output
    
    def annotate_dataset(self):
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
        """
        new_info = []
        for i, row in self.dataset.iterrows():
            acc = row[self.accession_col]
            pep = row[self.peptide_col]
            annotation = self.annotate_peptide(acc, pep)
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