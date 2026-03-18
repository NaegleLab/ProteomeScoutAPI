# **ProteomeScoutAPI**


A lightweight API to talk to ProteomeScout flatfiles in Python and annotate phosphoproteomic datasets with protein and site-specific information

See full [documentation](https://naeglelab.github.io/ProteomeScoutAPI/) for installation and use cases.

## About ProteomeScout
Version 3.0 - March 2026

Originally: By Alex Holehouse, Washington University in St. Louis Contact alex.holehouse@gmail.com or contribute at [https://github.com/alexholehouse](https://github.com/alexholehouse)

Current Version: Naegle Lab, University of Virginia [https://github.com/naegleLab](https://github.com/naegleLab)

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

## Usage

ProteomeScoutAPI allows for fast querying of protein-specific information, allowing for users to:
1. Explore properties of individual protein(s) of interest (domains, PTMs, structure, etc.)
2. Annotate a phosphoproteomic dataset with additional information about the proteins and sites in the dataset (e.g. domains, GO terms, etc.)
       
 




# Contributing code
The code here is incredibly simple, and the few methods presented give a good example of how one should parse the ProteomeScout records. If you're interested in adding the ability to parse out other information please go ahead and make a pull request. Tests would be appreciated too!
