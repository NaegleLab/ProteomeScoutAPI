#**ProteomeScoutAPI**


A lightweight API to talk to ProteomeScout flatfiles in Python. 

# about
Version 1.3 - February 2015 

By Alex Holehouse, Washington University in St. Louis ([Pappu Lab](http://pappulab.wustl.edu/))

Contact alex.holehouse@gmail.com or contribute at [https://github.com/alexholehouse](https://github.com/alexholehouse)

# overview
ProteomeScoutAPI is a Python module which can be used to connect to and parse ProteomeScout flatfiles. Specifically, the goal of this module is to allow anyone to interact with ProteomeScout data without the need to

1. Repeatedly query the ProteomeScout sever

2. Have any knowledge of SQL, or use an SQL-Python ORM

3. Facilitate rapid exploration of the ProteomeScout dataset

# usage

0. Import the `proteomeScoutAPI` module (this assumes the code is in your current directory)

        from proteomeScoutAPI import ProteomeScoutAPI

1. Create the API object by loading a flat file when you create the object. ProteomeScout flat files can be obtained from [https://proteomescout.wustl.edu/compendia](https://proteomescout.wustl.edu/compendia)


        PTM_API = ProteomeScoutAPI('<flat filename goes here>')
      
    This creates a PTM_API object, which provides a dynamic way to access     ProteomeScout data.

2. Query the API object using the functions. Available functions are

        PTM_API.get_phosphosites(ID)
        PTM_API.get_PTMs(ID)
        PTM_API.get_mutants(ID)
       
 
    ### `get_phososphosites(<ID>)` 
    
    This returns a list of 3-mer tuples, where each item in the list is a different phosphosite, and the three positions in each tuple associated with a phosphosite correspond to the position in the sequence, the residue being phosphorylated (Y/S/T) and the type. Note that position indexes from **1** (to maintain consistency with bioinformatics, not Python lists).
    
    If the ID passed cannot be found in the database, then a `-1` value is returned.
    
    If the ID passed as no phosphosites associated with it, then an empty list  (`[]`)is returned
    
    ### `get_PTMs(<ID>)` 
    
    This returns a list of 3-mer tuples, where each item in the list is a different post translational modification (PTM), including phosphosites. The three positions in each tuple associated with a PTM correspond to the position in the sequence, the residue being phosphorylated (Y/S/T) and the type. Note that, like with phosphosites, position indexes from **1** (to maintain consistency with bioinformatics, not Python lists).
    
    If the ID passed cannot be found in the database, then a `-1` value is returned.
    
    If the ID passed as no PTMs associated with it, then an empty list  (`[]`)is returned

    
    ### `get_mutations(<ID>)` 
    
    This returns a list of 4-mer tuples, where each item in the list is a different mutation, and the first three positions in each tuple associated with a mutations correspond to the position in sequence, the original residue residue type, and the new mutation residue type. The fourth position in the tuple is the label that is associated (e.g. pathogenic, non_pathogenic, or '')Note that, as always, positions index from **1** (to maintain consistency with bioinformatics, not Python lists).
    
    If the ID passed cannot be found in the database, then a `-1` value is returned.
    
    If the ID passed as no mutants associated with it, then an empty list  (`[]`)is returned



    The `uniqueKeys` is a list of the unique accession numbers to provide an easy way to loop over all the unique entries. 
    
    **NB:** IDs are redundant - i.e. several IDs point to the same protein record. This is is deliberate (and makes interfacing easy) but cycling over all the IDs in the API object would be incorrect and lead to double counting of some proteins in an inconsistent manner.


# examples

    from proteomeScoutAPI import ProteomeScoutAPI 
    
     
    ID = "O46631"
     
    PTM_API =         ProteomeScoutAPI("proteomescout_mammalia_20140831.tsv")
    
    print PTM_API.get_mutations(ID)
    
    PTM_API.get_PTMs(ID)
    
    # the following loop prints all the phosphosites in all the proteins
    for ID in PTM_API.uniqueKeys():
        print PTM_API.get_phosphosites(ID)


# contributing code
The code here is incredibly simple, and the few methods presented give a good example of how one should parse the ProteomeScout records. If you're interested in adding the ability to parse out other information please go ahead and make a pull request. Tests would be appreciated too!
