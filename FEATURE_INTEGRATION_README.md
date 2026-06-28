# Feature Integration Guide

This guide explains how to add a new dataset field so it propagates consistently across:

- Core API parsing and accessors in ProteomeScoutAPI
- Site-level PTM annotation in get_annotated_PTMs
- Species reference builds in SpeciesReferenceDataset
- ProteomicDataset annotation output

The current architecture is designed so most feature integration happens once in a shared site-annotation pipeline.

## Architecture Overview

The integration flow is:

1. Parse raw column data from each protein record
2. Add/extend accessor methods for feature retrieval
3. Add site-level mapping logic in one shared path
4. Reuse shared annotations in downstream builders

Key extension points:

- Header validation and optional column handling: proteomeScoutAPI/api.py (check file logic)
- Shared site annotation registry: ProteomeScoutAPI.POSITION_ANNOTATION_REGISTRY
- Shared position annotator: ProteomeScoutAPI._annotate_positions
- PTM table annotation: ProteomeScoutAPI.get_annotated_PTMs
- Species table construction: SpeciesReferenceDataset.build_protein_reference_dataset
- Peptide annotation output: ProteomicDataset.annotate_peptide

## Step-by-Step: Add a New Field

### 1. Decide column type

Choose one of these patterns:

- Per-protein value: one value per protein record
- Per-site feature ranges: feature name and start/stop boundaries mapped to PTM sites
- Per-site predictions keyed by residue position (like Spy-C)

### 2. Parse and expose accessors

In ProteomeScoutAPI:

- Read raw value from self.database[ID]
- Add a focused accessor that returns parsed values
- Return -1 when accession does not exist
- Handle missing/empty columns in a backward-compatible way

For optional columns:

- Keep them out of required header checks
- Include them in optional warnings
- Return empty/-1 when missing

### 3. Integrate into shared site annotation

If the field maps to residue positions, integrate via shared annotation logic.

For range-style features:

1. Ensure data can be normalized to tuples: e.g. (name, start, end, extra)
2. Add a rule in POSITION_ANNOTATION_REGISTRY with:
   - context_key
   - name_column
   - optional extra_column
   - optional boolean_column
   - include_boolean
3. Ensure _collect_site_annotation_context provides context_key data

For position-keyed predictions:

- Add per-position lookup in _annotate_positions
- Normalize output formatting in helper methods

This automatically propagates to:

- get_annotated_PTMs
- SpeciesReferenceDataset via annotated PTM reuse
- ProteomicDataset via _annotate_positions

### 4. Expose in downstream output schemas

If you want the new feature in exported dataframes:

- Add output columns in SpeciesReferenceDataset.OUTPUT_COLUMNS
- Add peptide annotation keys in ProteomicDataset.annotate_peptide output
- If the feature is sparse and species-specific, add those columns to OPTIONAL_ANNOTATION_COLUMNS so empty columns can be dropped from species exports

### 5. Keep behavior backward compatible

Compatibility requirements:

- Missing column in older datasets must not crash
- Empty values should produce empty strings/False/NaN where expected
- Table outputs must return empty DataFrame (not list) when output_format='table'

### 6. Add tests

Add/extend tests in tests/:

- Accessor parsing and lookup
- Training/missing edge cases for predictions
- Empty-record behavior
- Propagation into:
  - get_annotated_PTMs
  - SpeciesReferenceDataset output
  - ProteomicDataset annotate_peptide and annotate_dataset

Recommended command set: (assumes a conda environment like pscoutapi_dev has been setup allowing you to use developer mode)

conda run -n pscoutapi_dev python -m unittest discover -s tests -p 'test_spyc_predictions.py' -v
conda run -n pscoutapi_dev python -m unittest discover -s tests -p 'test_species_reference_dataset.py' -v
conda run -n pscoutapi_dev python -m unittest discover -s tests -p 'test_proteomic_dataset_spyc.py' -v

## Minimal Implementation Checklist

- Optional header handled safely
- Accessor method added
- Site annotation integrated through shared pipeline
- Downstream dataframe columns added (if needed)
- Sparse species columns dropped when fully empty (if needed)
- Unit tests added and passing

## Example Pattern (SpyC)

SpyC follows this pattern:

- Parse per-protein string entries with site:probability:class:confidence
- Expose protein-level and site-level accessors
- Normalize output values including Training case for missing probability with present confidence
- Emit site-level columns from shared annotator
- Reuse in both species reference and proteomic dataset annotation paths

This is the preferred pattern for future residue-level prediction fields.
