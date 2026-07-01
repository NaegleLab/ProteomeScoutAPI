Feature Integration
===================

The contributor guide for adding new fields that propagate through the API,
PTM annotations, species reference datasets, and proteomic dataset annotation
is available in the repository root:

`FEATURE_INTEGRATION_README.md <../../FEATURE_INTEGRATION_README.md>`_

This guide covers:

- Parser and accessor patterns for new columns
- Shared site-level annotation integration
- Propagation to `get_annotated_PTMs`, `SpeciesReferenceDataset`, and `ProteomicDataset`
- Backward compatibility for optional columns
- Recommended tests and validation workflow
