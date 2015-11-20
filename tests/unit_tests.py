import os.path
import time
import unittest
import random

from proteomeScoutAPI import ProteomeScoutAPI

class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.API = ProteomeScoutAPI('example.dat')

    def test_file_read(self):
        self.assertEqual(577, len(self.API.uniqueKeys))

    def test_get_GO(self):
        self.assertEqual(31, len(self.API.get_GO('Q9NYE7')))
        self.assertEqual('GO:0000790-nuclear chromatin', self.API.get_GO('Q9NYE7')[0])

    def test_get_PTMs(self):
        self.assertEqual(3, len(self.API.get_PTMs('Q9NYE7')))
        self.assertEqual('Acetylation', self.API.get_PTMs('Q9NYE7')[0][2])

    def test_get_pfam_domains(self):
        self.assertEqual(2, len(self.API.get_domains('Q9NYE7','pfam')))
        self.assertEqual('zf-H2C2_2', self.API.get_domains('Q9NYE7','pfam')[0][0])

    def test_get_nearby_PTMs(self):
        self.assertEqual(2, len(self.API.get_nearbyPTMs('Q9NYE7',5,10)))
        self.assertEqual(2, len(self.API.get_nearbyPTMs('Q9NYE7',5,13)))

    def test_get_mutations(self):
        self.assertEqual(1, len(self.API.get_mutations('Q9NYE7')))
        self.assertEqual('VAR_019971 (In dbSNP:rs3741665.)', self.API.get_mutations('Q9NYE7')[0][4])

    def test_get_phosphosites(self):
        self.assertEqual(32, len(self.API.get_phosphosites('117210')))
        self.assertEqual('880', self.API.get_phosphosites('117210')[10][0])
        
    def test_get_sequence(self):
        self.assertEqual('MAPTWGPGMV', self.API.get_sequence('117210')[0:10])

    def test_get_species(self):
        self.assertEqual('homo sapiens', self.API.get_species('117210'))

        
