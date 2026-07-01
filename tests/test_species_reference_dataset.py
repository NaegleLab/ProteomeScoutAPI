import os
import tempfile
import unittest

from proteomeScoutAPI import SpeciesReferenceDataset


class TestSpeciesReferenceDataset(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.builder = SpeciesReferenceDataset()

    def test_build_protein_reference_dataset_contains_expected_columns(self):
        reference_table = self.builder.build_protein_reference_dataset('O00629')

        self.assertFalse(reference_table.empty)
        self.assertEqual(self.builder.OUTPUT_COLUMNS, list(reference_table.columns))

        reference_row = reference_table.loc[reference_table['site'] == 'T24'].iloc[0]
        self.assertEqual('KPNA4', reference_row['gene_name'])
        self.assertEqual('Phosphothreonine', reference_row['ptm_type'])
        self.assertEqual(15, len(reference_row['oriented_peptide']))
        self.assertEqual('t', reference_row['oriented_peptide'][7])
        self.assertTrue(reference_row['in_interpro_domain'])
        self.assertIn('Importin-a_IBB', reference_row['interpro_domains'])
        self.assertIn('in_activation_loop', reference_table.columns)
        self.assertIn('activation_loop_qualities', reference_table.columns)

    def test_activation_loop_columns_are_consistent(self):
        reference_table = self.builder.build_protein_reference_dataset('O00629')

        expected_mask = reference_table['activation_loop_qualities'].fillna('').astype(str).str.len() > 0
        self.assertTrue((reference_table['in_activation_loop'] == expected_mask).all())

    def test_get_activation_loops_and_site_membership_from_new_format_string(self):
        accession = 'O00629'
        record = self.builder.database[accession]
        original_activation_loop = record.get('activation_loop', None)

        record['activation_loop'] = '20:30:high;200:210:medium'
        try:
            loops = self.builder.get_activation_loops(accession, output_format='list')
            self.assertEqual([('high', '20', '30'), ('medium', '200', '210')], loops)

            reference_table = self.builder.build_protein_reference_dataset(accession)
            reference_row = reference_table.loc[reference_table['site'] == 'T24'].iloc[0]
            self.assertTrue(reference_row['in_activation_loop'])
            self.assertEqual('high', reference_row['activation_loop_qualities'])
        finally:
            if original_activation_loop is None:
                record.pop('activation_loop', None)
            else:
                record['activation_loop'] = original_activation_loop

    def test_build_protein_reference_dataset_pads_protein_termini(self):
        reference_table = self.builder.build_protein_reference_dataset('O00401')

        reference_row = reference_table.loc[reference_table['site'] == 'S2'].iloc[0]
        self.assertEqual('______MsSVQQQPP', reference_row['oriented_peptide'])

    def test_build_species_reference_dataset_writes_csv(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            original_method = self.builder.return_species_nr_uniprot_ids
            self.builder.return_species_nr_uniprot_ids = lambda: (
                {'Homo sapiens': ['O00629', 'Q9UJX3']},
                {'Homo sapiens': True},
            )

            try:
                output_path = os.path.join(temp_dir, 'homo_sapiens_reference_dataset.csv')
                reference_table = self.builder.build_species_reference_dataset('Homo sapiens', output_file=output_path)
            finally:
                self.builder.return_species_nr_uniprot_ids = original_method

            self.assertTrue(os.path.exists(output_path))
            self.assertFalse(reference_table.empty)
            self.assertTrue((reference_table['species'] == 'Homo sapiens').all())
            self.assertIn('O00629', set(reference_table['uniprot_id']))
            self.assertIn('Q9UJX3', set(reference_table['uniprot_id']))

    def test_hidden_only_ptms_are_excluded_unless_requested(self):
        default_ptms = self.builder.get_PTMs('O43251')
        all_ptms = self.builder.get_PTMs('O43251', include_hidden=True)

        self.assertFalse(any(position == '380' and residue == 'Y' for position, residue, _ in default_ptms))
        self.assertTrue(any(position == '380' and residue == 'Y' for position, residue, _ in all_ptms))

        default_reference = self.builder.build_protein_reference_dataset('O43251')
        all_reference = self.builder.build_protein_reference_dataset('O43251', include_hidden=True)

        self.assertFalse((default_reference['site'] == 'Y380').any())
        self.assertTrue((all_reference['site'] == 'Y380').any())

    def test_species_reference_spyc_training_probability_mapping(self):
        accession = 'O00629'
        record = self.builder.database[accession]
        original_spyc_predictions = record.get('spyc_predictions', None)

        record['spyc_predictions'] = 'T24:NaN:NaN:High'
        try:
            reference_table = self.builder.build_protein_reference_dataset(accession)
            reference_row = reference_table.loc[reference_table['site'] == 'T24'].iloc[0]

            self.assertIn('SpY-C Prediction', reference_table.columns)
            self.assertIn('SpY-C Prediction Probability', reference_table.columns)
            self.assertEqual('High', reference_row['SpY-C Prediction'])
            self.assertEqual('Training', reference_row['SpY-C Prediction Probability'])
        finally:
            if original_spyc_predictions is None:
                record.pop('spyc_predictions', None)
            else:
                record['spyc_predictions'] = original_spyc_predictions

    def test_species_build_drops_empty_spyc_columns(self):
        original_method = self.builder.return_species_nr_uniprot_ids
        self.builder.return_species_nr_uniprot_ids = lambda: (
            {'Homo sapiens': ['O00629']},
            {'Homo sapiens': True},
        )

        try:
            reference_table = self.builder.build_species_reference_dataset('Homo sapiens')
            self.assertNotIn('SpY-C Prediction', reference_table.columns)
            self.assertNotIn('SpY-C Prediction Probability', reference_table.columns)
        finally:
            self.builder.return_species_nr_uniprot_ids = original_method