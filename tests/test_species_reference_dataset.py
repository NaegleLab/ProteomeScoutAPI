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