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


# Exception if file is bad
class BadProteomeScoutFile(Exception):
    pass


class ProteomeScoutAPI:


    def __init__(self,filename):
        """ 
        filename should be a ProteomeScout flatfile
        """        

        self.database={}
        self.uniqueKeys=[]

        # this will throw and exception if there's a problem with the file
        self.__checkFile(filename)
        
        self.__buildAPI(filename)

    def __checkFile(self, filename):
        """
        Internal function to check we (apparently) have a valid
        proteomeScout file
        """
        
        try:
            with open(filename, 'r') as f:
                first_line = f.readline()
                
            if not len(first_line.split("\t")) == 18:
                raise BadProteomeScoutFile("N/A")
            
                
        except:
            raise BadProteomeScoutFile("Invalid ProteomeScout flat file %s.\nFile is invalid or corrupted" % str(filename))
        

    def __buildAPI(self, datafile):

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
            for i in range(1,18):
                OBJ[headers[i]] = record[i].strip() if i < len(record) else ""

            # ALWAYS add the record via uniprot_id (or ProteomeScoutID if no uniprot_id)
            if uniprot_id in self.database:
                # if the newly found record has more PTMs associated with it than the 
                # original record then overwrite 
                if len(OBJ['modifications'].split(";")) > len(self.database[uniprot_id]['modifications'].split(";")):
                    self.database[uniprot_id] = OBJ
            else:
                self.database[uniprot_id] = OBJ

    def get_PTMs(self, ID):
        """
        Return all PTMs associated with the ID in question.

        POSTCONDITIONS:

        Returns a list of tuples of modifications
        [(position, residue, modification-type),...,]
        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications        

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        mods = record["modifications"]
        
        mods_raw=mods.split(";")
        mods_clean =[]
        for i in mods_raw:
            tmp = i.strip()
            tmp = tmp.split("-")
            
            # append a tuple of (position, residue, type)
            mods_clean.append((tmp[0][1:], tmp[0][0], "-".join(tmp[1:])))
        return mods_clean
    
    def get_structure(self, ID):
        """
        Return all structures associated with the ID in question.
       

        POSTCONDITIONS:
    
            Returns a list of tuples of structure
            if there is a problem with the start and end position, these will  be
            returned as -1
            [(domain_name, start_position, end_position),...,]
            
            Returns -1 if unable to find the ID
    
            Returns [] (empty list) if no structures  

        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
            
        structs = record["structure"]
               
        structs_raw=structs.split(";")
        structs_clean =[]
        for i in structs_raw:
            if i:
                tmp = i.strip()
                tmp = tmp.split(":")
                if len(tmp) >= 2:
                    structs_clean.append((tmp[0], tmp[1], tmp[2]))
        return structs_clean
    
    def get_macro_molecular(self, ID):
        """
        Return all macro-molecular structures associated with the ID in question.

        POSTCONDITIONS:
    
            Returns a list of tuples of macro-molecular structures
            [(structure_name, start_position, end_position),...,]
            
            Returns -1 if unable to find the ID
    
            Returns [] (empty list) if no macro-molecular structures  

        """
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
                    macro_mol_clean.append((name, start, stop))
                else:
                    print("ERROR: the macro-molecular structure did not match expected format %s"%(i))
        return macro_mol_clean
    
    def get_domains(self, ID, domain_type=None):
        """
        Return all domains associated with the ID in question.
        For interpro domains domain_type is 'interpro'
        For UniProt domains domain_type is 'uniprot'
        If domain_type is not specified, returns a dictionary with both.

        POSTCONDITIONS:

        Returns a list of tuples of domains if domain_type is specified (tuple length matched for easy processing, but only interpro_id is returned as the last for interpro):
        For interpro: [(domain_name, start_position, end_position, interpro_id),...,]
        For uniprot: [(domain_name, start_position, end_position, None),...,]
        
        Returns a dictionary with 'interpro' and 'uniprot' keys if domain_type is not specified:
        {'interpro': [...], 'uniprot': [...]}

        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no domains for the specified type

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        # If no domain_type specified, return both
        if domain_type is None:
            return {
                'interpro': self.get_domains(ID, 'interpro'),
                'uniprot': self.get_domains(ID, 'uniprot')
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
        return doms_clean
    
    def get_Scansite(self, ID):
        """
        DEPRECATED: Scansite predictions have been removed from the new data format.
        This method is kept for backward compatibility but will return an empty list.

        POSTCONDITIONS:

        Returns [] (empty list)

        """
        return []

    def get_Scansite_byPos(self, ID, res_pos):
        """
        DEPRECATED: Scansite predictions have been removed from the new data format.
        This method is kept for backward compatibility but will return an empty dictionary.

        POSTCONDITIONS:

        Returns {} (empty dictionary)

        """
        return {}

    def get_nearbyPTMs(self,ID,pos, window):
        """
        Return all PTMs within a specified window of a given position
        """
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
        return nearby_mods

    def get_species(self,ID):
        """
        Return the species associated with the ID
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return record["species"]
    
    def get_sequence(self, ID):
        """
        Return the sequence associated with the ID
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return record["sequence"]

    def get_acc_gene(self, ID):
        """
        Return the accession number and gene name associated with the ID
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        return (record["accession"], record["gene"])
         

    
    def get_phosphosites(self,ID):
        """
        Return all phosphosites (S/T/Y phosphorylation) associated with the ID
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        
        mods = self.get_PTMs(ID)
        phos_sites = []
        for mod in mods:
            residue = mod[1]
            mod_type = mod[2]
            if mod_type in ['Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine']:
                phos_sites.append(mod)
        return phos_sites

    def get_evidence(self, ID):
        """
        Return evidence scores associated with PTMs for the ID
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

        POSTCONDITIONS:
        
        Returns a list of tuples (GO_term, type) where type is 'F' (Molecular Function), 
        'P' (Biological Process), or 'C' (Cellular Component)

        Returns a -1 if unable to find the ID

        Returns a [] (empty list) if no GO terms

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
        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        accessions = record["accession"].split(";")
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

        POSTCONDITIONS:
        Returns a tuple (species_dict, species_reference_bool)
        species_dict: {species_name: [list of uniprot IDs]}
        species_reference_bool: {species_name: is_reference_species (True/False)} - If False, it means that only nonredundant uniprot records are included in the list if they have PTMs. 
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