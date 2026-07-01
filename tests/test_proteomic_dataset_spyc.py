import unittest
from unittest.mock import patch

import pandas as pd

from proteomeScoutAPI import ProteomicDataset


class TestProteomicDatasetSpyC(unittest.TestCase):

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
