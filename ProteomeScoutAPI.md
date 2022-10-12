# Table of Contents

* [proteomeScoutAPI](#proteomeScoutAPI)
  * [ProteomeScoutAPI](#proteomeScoutAPI.ProteomeScoutAPI)
    * [\_\_init\_\_](#proteomeScoutAPI.ProteomeScoutAPI.__init__)
    * [get\_PTMs](#proteomeScoutAPI.ProteomeScoutAPI.get_PTMs)
    * [get\_structure](#proteomeScoutAPI.ProteomeScoutAPI.get_structure)
    * [get\_domains](#proteomeScoutAPI.ProteomeScoutAPI.get_domains)
    * [get\_domains\_harmonized](#proteomeScoutAPI.ProteomeScoutAPI.get_domains_harmonized)
    * [get\_Scansite](#proteomeScoutAPI.ProteomeScoutAPI.get_Scansite)
    * [get\_Scansite\_byPos](#proteomeScoutAPI.ProteomeScoutAPI.get_Scansite_byPos)
    * [get\_nearbyPTMs](#proteomeScoutAPI.ProteomeScoutAPI.get_nearbyPTMs)
    * [get\_species](#proteomeScoutAPI.ProteomeScoutAPI.get_species)
    * [get\_sequence](#proteomeScoutAPI.ProteomeScoutAPI.get_sequence)
    * [get\_acc\_gene](#proteomeScoutAPI.ProteomeScoutAPI.get_acc_gene)
    * [get\_phosphosites](#proteomeScoutAPI.ProteomeScoutAPI.get_phosphosites)
    * [get\_evidence](#proteomeScoutAPI.ProteomeScoutAPI.get_evidence)
    * [get\_mutations](#proteomeScoutAPI.ProteomeScoutAPI.get_mutations)
    * [get\_GO](#proteomeScoutAPI.ProteomeScoutAPI.get_GO)
    * [get\_kinaseLoops](#proteomeScoutAPI.ProteomeScoutAPI.get_kinaseLoops)
    * [get\_accessions](#proteomeScoutAPI.ProteomeScoutAPI.get_accessions)
    * [get\_PTMs\_withEvidenceThreshold](#proteomeScoutAPI.ProteomeScoutAPI.get_PTMs_withEvidenceThreshold)
    * [get\_PTMs\_withEvidence](#proteomeScoutAPI.ProteomeScoutAPI.get_PTMs_withEvidence)

<a id="proteomeScoutAPI"></a>

# proteomeScoutAPI

<a id="proteomeScoutAPI.ProteomeScoutAPI"></a>

## ProteomeScoutAPI Objects

```python
class ProteomeScoutAPI()
```

<a id="proteomeScoutAPI.ProteomeScoutAPI.__init__"></a>

#### \_\_init\_\_

```python
def __init__(filename)
```

filename should be a ProteomeScout flatfile

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_PTMs"></a>

#### get\_PTMs

```python
def get_PTMs(ID)
```

Return all PTMs associated with the ID in question.

POSTCONDITIONS:

Returns a list of tuples of modifications
[(position, residue, modification-type),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no modifications

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_structure"></a>

#### get\_structure

```python
def get_structure(ID)
```

Return all structures associated with the ID in question.


POSTCONDITIONS:

    Returns a list of tuples of structure
    if there is a problem with the start and end position, these will  be
    returned as -1
    [(domain_name, start_position, end_position),...,]

    Returns -1 if unable to find the ID

    Returns [] (empty list) if no structures

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_domains"></a>

#### get\_domains

```python
def get_domains(ID, domain_type)
```

Return all domains associated with the ID in question.
For pfam domains domain_type is 'pfam'
For UniProt domains domain_type is 'uniprot'

POSTCONDITIONS:

Returns a list of tuples of domains 
if there is a problem with the start and end position, these will be
returned as -1
[(domain_name, start_position, end_position),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no modifications

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_domains_harmonized"></a>

#### get\_domains\_harmonized

```python
def get_domains_harmonized(ID)
```

This will harmonize pfam and uniprot domains into one tuple output
with the following rules:
    1. Use the pfam domain name (since it does not append numbers if more than one domain of that type)
    2. Use the Uniprot domain boundaries (since they are typically more expansive)
    3. If uniprot does not have a pfam domain, use the pfam domain as it is
POSTCONDITIONS:

Returns a list of tuples of domains 
if there is a problem with the start and end position, these will be
returned as -1
[(domain_name, start_position, end_position),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no domains

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_Scansite"></a>

#### get\_Scansite

```python
def get_Scansite(ID)
```

Return all Scansite annotations that exist for a protein record (ID)


POSTCONDITIONS:

Returns a dict of tuples of scansite predictions. dictionary keys are the type of Scansite prediction (bind, kinase)
Tuples in each are 
[(residuePosition, kinase/bind name, score),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no scansite predictions

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_Scansite_byPos"></a>

#### get\_Scansite\_byPos

```python
def get_Scansite_byPos(ID, res_pos)
```

Return all Scansite annotations that exist at a specific residue position, given as the AminoAcidPos, e.g. T183

POSTCONDITIONS:

Returns a list of tuples of modifications
[(position, residue, modification-type),...,]

Returns -1 if unable to find the ID

Returns empty dictionary if no Scansite. Returns emtpy lists in Dictionary items
if no Scansite predictions found at that site

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_nearbyPTMs"></a>

#### get\_nearbyPTMs

```python
def get_nearbyPTMs(ID, pos, window)
```

Return all PTMs associated with the ID in question that reside within
+/-window, relative to the designated position (pos)

POSTCONDITIONS:

Returns a list of tuples of modifications
[(position, residue, modification-type),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no modifications

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_species"></a>

#### get\_species

```python
def get_species(ID)
```

Return the species associated with the ID in question.
POSTCONDITIONS:

Returns a string of the species name

Returns '-1' if unable to find the ID

Returns '' (empty list) if no species

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_sequence"></a>

#### get\_sequence

```python
def get_sequence(ID)
```

Return the sequence associated with the ID in question.
POSTCONDITIONS:

Returns a string of the sequence

Returns '-1' if unable to find the ID

Returns '' (empty list) if no sequence

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_acc_gene"></a>

#### get\_acc\_gene

```python
def get_acc_gene(ID)
```

Return the gene name  associated with the ID in question.
POSTCONDITIONS:

Returns a string of the gene name 

Returns '-1' if unable to find the ID

Returns '' (empty list) if no gene name

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_phosphosites"></a>

#### get\_phosphosites

```python
def get_phosphosites(ID)
```

Return all phosphosites associated with the ID in question.

POSTCONDITIONS:

Returns a list of tuples of phosphosites
[(position, residue, phosphosite-type),...,]

Returns -1 if unable to find the ID

Returns [] (empty list) if no modifications

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_evidence"></a>

#### get\_evidence

```python
def get_evidence(ID)
```

Return all evidence associated with modifications for the ID in question.

POSTCONDITIONS:

Returns a list of tuples of evidence
[(numEvidences, [ArrayOfEvidences], (numEvidences ScecondSite, [ArrayOfEvidences], etc.]

Returns -1 if unable to find the ID

Returns [] (empty list) if no modifications/evidences

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_mutations"></a>

#### get\_mutations

```python
def get_mutations(ID)
```

Return all mutations associated with the ID in question.
mutations = PTM_API.get_mutations(ID)

POSTCONDITIONS:

Returns a list of tuples of mutations 
[(original residue, position, new residue, annotation),...,]

Returns -1 if unable to find the ID
Returns -2 if the number of mutations and annotations do not match

Returns [] (empty list) if no mutations

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_GO"></a>

#### get\_GO

```python
def get_GO(ID)
```

Return all GO terms associated with the ID in question

POSTCONDITIONS:

Returns a list of GO Terms

Returns a -1 if unable to find the ID

Returns a [] (empty list) if no GO terms

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_kinaseLoops"></a>

#### get\_kinaseLoops

```python
def get_kinaseLoops(ID)
```

Return kinase activation loop with the ID in question


POSTCONDITIONS:

Returns a tuple of (loop/predicted, start, stop)

Returns a -1 if unable to find the ID

Returns a [] (empty list) if no kinase activation loops

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_accessions"></a>

#### get\_accessions

```python
def get_accessions(ID)
```

Return a list of accessions associated with the protein

POSTCONDITIONS:

Returns a list of the protein's accessions

Returns a -1 if unable to find the ID

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_PTMs_withEvidenceThreshold"></a>

#### get\_PTMs\_withEvidenceThreshold

```python
def get_PTMs_withEvidenceThreshold(ID, evidenceThreshold)
```

Return all PTMs associated with the ID in question that have at least evidenceThreshold or more pieces
of experimental or database evidence

POSTCONDITIONS:

Returns a list of tuples of modifications
[(position, residue, modification-type),...,]

Returns -1 if unable to find the ID
Returns -2 if the modifications have a mismatched evidence array

Returns [] (empty list) if no modifications fit that threshold

<a id="proteomeScoutAPI.ProteomeScoutAPI.get_PTMs_withEvidence"></a>

#### get\_PTMs\_withEvidence

```python
def get_PTMs_withEvidence(ID)
```

Return PTMs as a list of dictionaries with both the modification information and the evidence associated with that.

POSTCONDITIONS:

Returns a list of dictionaries, keys are 'mod' and 'evidence'.
dictionary[0]['mod'] is an array with (position, residue, modification-type) of the first PTM in the ID record
dictionary[0]['evidence'] is a list of experiment IDs for that PTM

Returns -1 if unable to find ID
Returns -2 if the modifications have a mimatched evidence array

Returns [] (empty list) if no modifications

