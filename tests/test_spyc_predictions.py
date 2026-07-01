import unittest

from proteomeScoutAPI import SpeciesReferenceDataset


class TestSpyCPredictions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.api = SpeciesReferenceDataset()

    def test_get_spyc_predictions_and_site_lookup(self):
        accession = 'O43251'
        record = self.api.database[accession]
        original_spyc_predictions = record.get('spyc_predictions', None)

        record['spyc_predictions'] = 'Y380:0.82:positive:High;S381:NaN:NaN:Medium;T2:NaN:NaN:NaN'
        try:
            predictions = self.api.get_spyc_predictions(accession)
            self.assertIsInstance(predictions, dict)
            self.assertEqual(['0.82', 'positive', 'High'], predictions[380])
            self.assertEqual(['NaN', 'NaN', 'Medium'], predictions[381])

            self.assertEqual(('0.82', 'positive', 'High'), self.api.get_spyc_predictions_byPos(accession, 380))
            self.assertEqual(('NaN', 'NaN', 'Medium'), self.api.get_spyc_predictions_byPos(accession, 381))
            self.assertEqual(-1, self.api.get_spyc_predictions_byPos(accession, 999999))
        finally:
            if original_spyc_predictions is None:
                record.pop('spyc_predictions', None)
            else:
                record['spyc_predictions'] = original_spyc_predictions

    def test_spyc_backward_compatible_when_missing(self):
        accession = 'O43251'
        record = self.api.database[accession]
        had_key = 'spyc_predictions' in record
        original_spyc_predictions = record.get('spyc_predictions', None)

        if had_key:
            record.pop('spyc_predictions', None)
        try:
            self.assertEqual(-1, self.api.get_spyc_predictions(accession))
            self.assertEqual(-1, self.api.get_spyc_predictions_byPos(accession, 380))

            annotated = self.api.get_annotated_PTMs(accession, include_hidden=True)
            self.assertIn('SpY-C Prediction', annotated.columns)
            self.assertIn('SpY-C Prediction Probability', annotated.columns)
            self.assertTrue((annotated['SpY-C Prediction'].fillna('') == '').all())
            self.assertTrue((annotated['SpY-C Prediction Probability'].fillna('') == '').all())
        finally:
            if had_key:
                record['spyc_predictions'] = original_spyc_predictions

    def test_training_label_in_annotated_ptms(self):
        accession = 'O43251'
        record = self.api.database[accession]
        original_spyc_predictions = record.get('spyc_predictions', None)

        record['spyc_predictions'] = 'Y380:NaN:NaN:Medium'
        try:
            annotated = self.api.get_annotated_PTMs(accession, include_hidden=True)
            site_row = annotated.loc[annotated['Position'].astype(str) == '380']
            self.assertFalse(site_row.empty)

            row = site_row.iloc[0]
            self.assertEqual('Medium', row['SpY-C Prediction'])
            self.assertEqual('Training', row['SpY-C Prediction Probability'])
        finally:
            if original_spyc_predictions is None:
                record.pop('spyc_predictions', None)
            else:
                record['spyc_predictions'] = original_spyc_predictions

    def test_get_ptms_table_empty_returns_dataframe(self):
        accession = 'O43251'
        record = self.api.database[accession]
        original_modifications = record.get('modifications', None)

        record['modifications'] = ''
        try:
            ptms_table = self.api.get_PTMs(accession, output_format='table', include_hidden=True)
            self.assertTrue(hasattr(ptms_table, 'empty'))
            self.assertTrue(ptms_table.empty)

            annotated = self.api.get_annotated_PTMs(accession, include_hidden=True)
            self.assertTrue(hasattr(annotated, 'empty'))
            self.assertTrue(annotated.empty)
        finally:
            if original_modifications is None:
                record.pop('modifications', None)
            else:
                record['modifications'] = original_modifications
