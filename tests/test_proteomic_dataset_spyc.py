import unittest
from unittest.mock import patch

import pandas as pd

from proteomeScoutAPI import ProteomicDataset, helpers


class TestProteomicDatasetSpyC(unittest.TestCase):

    def test_lowercase_lysine_is_detected_as_modified_site(self):
        self.assertEqual([3], helpers.find_mod('AAAkAAAK'))

        aligned_peptides, positions, residues = helpers.returnOrientedPhosphoPeptide('AAAKAAAK', 'AAAkAAAK', flank=1)

        self.assertEqual([3], positions)
        self.assertEqual(['K'], residues)
        self.assertEqual(1, len(aligned_peptides))
        self.assertIn('k', aligned_peptides[0])

    def test_annotate_peptide_includes_spyc_columns(self):
        dataset = pd.DataFrame({'acc': ['O43251'], 'pep': ['YTEST']})
        annotator = ProteomicDataset(dataset, accession_col='acc', peptide_col='pep', find_site=True, update=False)

        record = annotator.database['O43251']
        original_spyc_predictions = record.get('spyc_predictions', None)
        record['spyc_predictions'] = 'Y380:NaN:NaN:Medium'

        try:
            with patch('proteomeScoutAPI.api.helpers.returnOrientedPhosphoPeptide', return_value=(['PEPTIDE'], [380], ['Y'])):
                annotation = annotator.annotate_peptide('O43251', 'YTEST', include_hidden=True)

            self.assertEqual('Medium', annotation['SpY-C Prediction'])
            self.assertEqual('Training', annotation['SpY-C Prediction Probability'])
        finally:
            if original_spyc_predictions is None:
                record.pop('spyc_predictions', None)
            else:
                record['spyc_predictions'] = original_spyc_predictions

    def test_annotate_dataset_appends_spyc_columns(self):
        dataset = pd.DataFrame({'acc': ['O43251'], 'pep': ['YTEST']})
        annotator = ProteomicDataset(dataset, accession_col='acc', peptide_col='pep', find_site=True, update=False)

        record = annotator.database['O43251']
        original_spyc_predictions = record.get('spyc_predictions', None)
        record['spyc_predictions'] = 'Y380:0.83:positive:High'

        try:
            with patch('proteomeScoutAPI.api.helpers.returnOrientedPhosphoPeptide', return_value=(['PEPTIDE'], [380], ['Y'])):
                annotator.annotate_dataset(include_hidden=True)

            self.assertIn('SpY-C Prediction', annotator.dataset.columns)
            self.assertIn('SpY-C Prediction Probability', annotator.dataset.columns)
            self.assertEqual('High', annotator.dataset.loc[0, 'SpY-C Prediction'])
            self.assertEqual('0.83', annotator.dataset.loc[0, 'SpY-C Prediction Probability'])
        finally:
            if original_spyc_predictions is None:
                record.pop('spyc_predictions', None)
            else:
                record['spyc_predictions'] = original_spyc_predictions

    def test_documented_modification_details_commas_multiple_types_same_site(self):
        dataset = pd.DataFrame({'acc': ['O43251'], 'pep': ['YTEST']})
        annotator = ProteomicDataset(dataset, accession_col='acc', peptide_col='pep', find_site=True, update=False)

        record = annotator.database['O43251']
        original_modifications = record.get('modifications', None)
        original_evidence = record.get('evidence', None)
        record['modifications'] = 'Y380-Phosphotyrosine;Y380-Acetylation'
        record['evidence'] = ''

        try:
            with patch('proteomeScoutAPI.api.helpers.returnOrientedPhosphoPeptide', return_value=(['PEPTIDE'], [380], ['Y'])):
                annotation = annotator.annotate_peptide('O43251', 'YTEST', include_hidden=False)

            self.assertEqual('1', annotation['documented_modification'])
            self.assertIn('Phosphotyrosine, Acetylation', annotation['documented_modification_details'])
            self.assertNotIn('Phosphotyrosine;Acetylation', annotation['documented_modification_details'])
        finally:
            if original_modifications is None:
                record.pop('modifications', None)
            else:
                record['modifications'] = original_modifications

            if original_evidence is None:
                record.pop('evidence', None)
            else:
                record['evidence'] = original_evidence
