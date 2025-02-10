import unittest
from unittest.mock import MagicMock, mock_open, patch
import numpy as np
from osmofold import osmofold_local
import mdtraj

class TestGetMaxSasaList(unittest.TestCase):
    
    def test_return_type(self):
        """Test that the function returns a dictionary."""
        result = osmofold_local.get_max_sasa_list()
        self.assertIsInstance(result, dict)
    
    def test_keys_exist(self):
        """Test that the dictionary contains both 'backbone' and 'sidechain' keys."""
        result = osmofold_local.get_max_sasa_list()
        self.assertIn("backbone", result)
        self.assertIn("sidechain", result)
    
    def test_inner_structure(self):
        """Test that both 'backbone' and 'sidechain' contain dictionaries."""
        result = osmofold_local.get_max_sasa_list()
        self.assertIsInstance(result["backbone"], dict)
        self.assertIsInstance(result["sidechain"], dict)
    
    def test_expected_values(self):
        """Test that specific values are correct."""
        result = osmofold_local.get_max_sasa_list()
        self.assertAlmostEqual(result["backbone"]["A"], 27.85)
        self.assertAlmostEqual(result["sidechain"]["R"], 171.10)
    
    def test_all_amino_acids_present(self):
        """Test that all 20 amino acids are present in both backbone and sidechain."""
        expected_amino_acids = set("ARNDCEQGHILKMFPSTWYV")
        result = osmofold_local.get_max_sasa_list()
        self.assertEqual(set(result["backbone"].keys()), expected_amino_acids)
        self.assertEqual(set(result["sidechain"].keys()), expected_amino_acids)

class TestAminoToEnergy(unittest.TestCase):
    def test_return_type(self):
        self.assertIsInstance(osmofold_local.amino_to_energy("A", "trehalose_hong"), float)
    
    def test_known_values(self):
        self.assertAlmostEqual(osmofold_local.amino_to_energy("A", "trehalose_hong"), 59.3)
        self.assertAlmostEqual(osmofold_local.amino_to_energy("R", "tmao"), -109.27)
    
    def test_back_value(self):
        self.assertEqual(osmofold_local.amino_to_energy("A", "trehalose_hongBack"), 35)
        self.assertEqual(osmofold_local.amino_to_energy("R", "tmaoBack"), 90)
    
    def test_invalid_amino_acid(self):
        self.assertIsNone(osmofold_local.amino_to_energy("X", "trehalose_hong"))
    
    def test_invalid_cosolute(self):
        self.assertIsNone(osmofold_local.amino_to_energy("A", "unknown"))

class TestExtractSequence(unittest.TestCase):

    def test_return_type(self):
        """Test that the function returns a string."""
        mock_data = """ATOM      1  N   MET A   1      11.104  13.207   9.657  1.00 20.00           N  \n"""
        with patch("builtins.open", mock_open(read_data=mock_data)):
            self.assertIsInstance(osmofold_local.extract_sequence("example.pdb"), str)
    
    def test_empty_file(self):
        """Test that an empty file returns an empty string."""
        with open("empty.pdb", "w") as f:
            pass
        self.assertEqual(osmofold_local.extract_sequence("empty.pdb"), "")
    
    def test_known_sequence(self):
        """Test extraction with a known PDB input."""
        pdb_content = """\nATOM      1  N   MET A   1      11.104  13.207   9.657  1.00 20.00           N  \nATOM      2  CA  MET A   1      12.560  13.207   9.657  1.00 20.00           C  \nATOM      3  N   ALA A   2      14.104  13.207   9.657  1.00 20.00           N  \n"""
        with open("test.pdb", "w") as f:
            f.write(pdb_content)
        self.assertEqual(osmofold_local.extract_sequence("test.pdb"), "MA")
    
    def test_no_atoms(self):
        """Test that a PDB file without ATOM/HETATM lines returns an empty string."""
        pdb_content = """\nHEADER    TEST PDB\nTITLE     EXAMPLE WITHOUT ATOMS\nEND\n"""
        with open("no_atoms.pdb", "w") as f:
            f.write(pdb_content)
        self.assertEqual(osmofold_local.extract_sequence("no_atoms.pdb"), "")

class TestThreetoOne(unittest.TestCase):

    def test_three_to_one_valid_residues(self):
        """
        Test the three_to_one function for valid three-letter residue codes.
        """
        test_cases = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }

        for three_letter, one_letter in test_cases.items():
            with self.subTest(residue=three_letter):
                result = osmofold_local.three_to_one(three_letter)
                self.assertEqual(result, one_letter)

    def test_three_to_one_invalid_residues(self):
        """
        Test the three_to_one function for invalid or unknown three-letter residue codes.
        """
        invalid_residues = ['XYZ', 'ZZZ', '123', '', 'AAR']
        for residue in invalid_residues:
            with self.subTest(residue=residue):
                result = osmofold_local.three_to_one(residue)
                self.assertEqual(result, 'X')  # 'X' is returned for unknown residues

    def test_three_to_one_case_insensitivity(self):
        """
        Test the three_to_one function to ensure it handles case insensitivity.
        """
        result = osmofold_local.three_to_one('ala'.upper())  # Convert to uppercase
        self.assertEqual(result, 'A')

        result = osmofold_local.three_to_one('arg'.upper())  # Convert to uppercase
        self.assertEqual(result, 'R')

class TestExtractSequencesByChains(unittest.TestCase):

    @patch("osmofold.osmofold_local.three_to_one", side_effect=lambda x: x[0])  # Mock three_to_one to return first letter
    def test_single_chain(self, mock_three_to_one):
        pdb_data = """\
ATOM      1  N   MET A   1      
ATOM      2  CA  MET A   1      
ATOM      3  C   MET A   1      
ATOM      4  N   GLY A   2      
ATOM      5  CA  GLY A   2      
ATOM      6  C   GLY A   2      
"""
        with open("test.pdb", "w") as f:
            f.write(pdb_data)

        result = osmofold_local.extract_sequences_by_chains("test.pdb")
        self.assertEqual(result, ['MG'])

    @patch("osmofold.osmofold_local.three_to_one", side_effect=lambda x: x[0])
    def test_multiple_chains(self, mock_three_to_one):
        pdb_data = """\
ATOM      1  N   ALA A   1      
ATOM      2  CA  ALA A   1      
ATOM      3  C   ALA A   1      
ATOM      4  N   LYS B   1      
ATOM      5  CA  LYS B   1      
ATOM      6  C   LYS B   1      
ATOM      7  N   GLY B   2      
ATOM      8  CA  GLY B   2      
ATOM      9  C   GLY B   2      
"""
        with open("test.pdb", "w") as f:
            f.write(pdb_data)

        result = osmofold_local.extract_sequences_by_chains("test.pdb")
        self.assertEqual(result, ['A', 'LG'])

    @patch("osmofold.osmofold_local.three_to_one", side_effect=lambda x: x[0])
    def test_duplicate_residues_ignored(self, mock_three_to_one):
        pdb_data = """\
ATOM      1  N   SER A   1      
ATOM      2  CA  SER A   1      
ATOM      3  C   SER A   1      
ATOM      4  O   SER A   1      
ATOM      5  N   GLY A   2      
ATOM      6  CA  GLY A   2      
"""
        with open("test.pdb", "w") as f:
            f.write(pdb_data)

        result = osmofold_local.extract_sequences_by_chains("test.pdb")
        self.assertEqual(result, ['SG'])

    @patch("osmofold.osmofold_local.three_to_one", side_effect=lambda x: 'X' if x not in ['ALA', 'GLY', 'LYS'] else x[0])
    def test_invalid_residue_name(self, mock_three_to_one):
        pdb_data = """\
ATOM      1  N   XYZ A   1      
ATOM      2  CA  XYZ A   1      
ATOM      3  C   XYZ A   1      
ATOM      4  N   ALA A   2      
ATOM      5  CA  ALA A   2      
"""
        with open("test.pdb", "w") as f:
            f.write(pdb_data)

        result = osmofold_local.extract_sequences_by_chains("test.pdb")
        self.assertEqual(result, ['XA'])

class TestGetTFE(unittest.TestCase):

    @patch('osmofold.osmofold_local.amino_to_energy')
    def test_get_tfe_with_custom_tfe(self, mock_amino_to_energy):
        custom_tfe = {
            "backbone": {'A': 1.0, 'F': 2.0, 'L': 1.5},
            "sidechain": {'A': 0.5, 'F': 1.0, 'L': 0.8}
        }

        # Do not mock amino_to_energy here since custom_tfe should be used directly
        # Call the function with custom TFE
        seq = 'AF'
        osmo = 'osmolyte'
        backbone_tfe, sidechain_tfe = osmofold_local.get_tfe(seq, osmo, custom_tfe)

        # Assertions
        self.assertEqual(backbone_tfe, [1.0, 2.0])
        self.assertEqual(sidechain_tfe, [0.5, 1.0])

    @patch('osmofold.osmofold_local.amino_to_energy')
    def test_get_tfe_without_custom_tfe(self, mock_amino_to_energy):
        # Mock amino_to_energy to return fixed values based on the osmolyte
        mock_amino_to_energy.side_effect = lambda aa, osmo: 0.1 if "Back" in osmo else 0.2
        
        seq = 'AF'
        osmo = 'osmolyte'
        
        # Call the function
        backbone_tfe, sidechain_tfe = osmofold_local.get_tfe(seq, osmo)
        
        # Assertions (values based on the mock behavior)
        self.assertEqual(backbone_tfe, [0.1, 0.1])  # Backbone TFE will be 0.1 for both A and F
        self.assertEqual(sidechain_tfe, [0.2, 0.2])  # Sidechain TFE will be 0.2 for both A and F

    @patch('osmofold.osmofold_local.amino_to_energy')
    def test_get_tfe_with_invalid_custom_tfe(self, mock_amino_to_energy):
        # Custom TFE with missing amino acid entries
        custom_tfe = {
            "backbone": {'A': 1.0, 'F': 2.0},
            "sidechain": {'A': 0.5, 'F': 1.0}
        }

        # Mock behavior to simulate missing amino acids
        mock_amino_to_energy.side_effect = lambda aa, osmo: 0.0

        seq = 'AF'
        osmo = 'osmolyte'
        
        # Call the function
        backbone_tfe, sidechain_tfe = osmofold_local.get_tfe(seq, osmo, custom_tfe)
        
        # Assertions (since custom TFE is missing some amino acids, we fall back to default)
        self.assertEqual(backbone_tfe, [1.0, 2.0])
        self.assertEqual(sidechain_tfe, [0.5, 1.0])

    @patch('osmofold.osmofold_local.amino_to_energy')
    def test_get_tfe_with_empty_sequence(self, mock_amino_to_energy):
        # Mock behavior for an empty sequence
        mock_amino_to_energy.side_effect = lambda aa, osmo: 0.0
        
        seq = ''
        osmo = 'osmolyte'
        
        # Call the function
        backbone_tfe, sidechain_tfe = osmofold_local.get_tfe(seq, osmo)
        
        # Assertions (no amino acids, so TFE lists should be empty)
        self.assertEqual(backbone_tfe, [])
        self.assertEqual(sidechain_tfe, [])

    @patch('osmofold.osmofold_local.amino_to_energy')
    def test_get_tfe_with_unexpected_aa(self, mock_amino_to_energy):
        # Custom TFE with unexpected amino acid
        custom_tfe = {
            "backbone": {'A': 1.0, 'F': 2.0, 'L': 1.5},
            "sidechain": {'A': 0.5, 'F': 1.0, 'L': 0.8}
        }

        # Mock amino_to_energy to always return 0 for non-existent amino acids
        mock_amino_to_energy.side_effect = lambda aa, osmo: 0.0

        seq = 'XYZ'  # Invalid amino acids
        osmo = 'osmolyte'
        
        # Call the function
        backbone_tfe, sidechain_tfe = osmofold_local.get_tfe(seq, osmo, custom_tfe)
        
        # Assertions (since custom TFE doesn't have X, Y, or Z, both lists will be [0, 0, 0])
        self.assertEqual(backbone_tfe, [0, 0, 0])
        self.assertEqual(sidechain_tfe, [0, 0, 0])

class TestSasaToRasa(unittest.TestCase):
    @patch("osmofold.osmofold_local.get_max_sasa_list")
    def test_sasa_to_rasa(self, mock_get_max_sasa_list):
        # Mock the max SASA values
        mock_get_max_sasa_list.return_value = {
            "backbone": {"A": 10.0, "C": 12.0},
            "sidechain": {"A": 20.0, "C": 25.0}
        }
        
        seq = "AC"
        backbone_sasa = [5.0, 6.0]
        sidechain_sasa = [10.0, 12.5]
        
        expected_backbone_rasa = [5.0 / 10.0, 6.0 / 12.0]  # [0.5, 0.5]
        expected_sidechain_rasa = [10.0 / 20.0, 12.5 / 25.0]  # [0.5, 0.5]
        
        backbone_rasa, sidechain_rasa = osmofold_local.sasa_to_rasa(seq, backbone_sasa, sidechain_sasa)
        
        self.assertEqual(backbone_rasa, expected_backbone_rasa)
        self.assertEqual(sidechain_rasa, expected_sidechain_rasa)

    def test_sasa_to_rasa_empty_sequence(self):

        seq = ""
        backbone_sasa = []
        sidechain_sasa = []
        
        backbone_rasa, sidechain_rasa = osmofold_local.sasa_to_rasa(seq, backbone_sasa, sidechain_sasa)
        
        self.assertEqual(backbone_rasa, [])
        self.assertEqual(sidechain_rasa, [])

class TestGetPDBInfo(unittest.TestCase):
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_valid_pdb(self, MockSSTrajectory):
        mock_pdb = MagicMock()
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = "ACDE"
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_all_SASA.side_effect = [[[1.1, 2.2, 3.3, 4.4]], [[0.5, 1.5, 2.5, 3.5]]]
        MockSSTrajectory.return_value = mock_pdb
        
        result = osmofold_local.get_pdb_info("dummy.pdb")
        self.assertEqual(result, ("ACDE", [1.1, 2.2, 3.3, 4.4], [0.5, 1.5, 2.5, 3.5]))
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_empty_pdb(self, MockSSTrajectory):
        mock_pdb = MagicMock()
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = ""
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_all_SASA.side_effect = [[], []]
        MockSSTrajectory.return_value = mock_pdb
        
        with self.assertRaises(Exception) as context:
            osmofold_local.get_pdb_info("invalid.pdb")
        self.assertEqual(str(context.exception), "list index out of range")

    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_invalid_pdb(self, MockSSTrajectory):
        MockSSTrajectory.side_effect = Exception("Invalid PDB format")
        
        with self.assertRaises(Exception) as context:
            osmofold_local.get_pdb_info("invalid.pdb")
        self.assertEqual(str(context.exception), "Invalid PDB format")

class TestGetChainInfo(unittest.TestCase):
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_valid_pdb_single_chain(self, MockSSTrajectory):
        mock_pdb = MagicMock()
        
        # Create a separate MagicMock for proteins
        mock_proteins = MagicMock()
        mock_proteins.get_amino_acid_sequence.return_value = "ACDE"
        mock_proteins.get_all_SASA.side_effect = [[[1.1, 2.2, 3.3, 4.4]], [[0.5, 1.5, 2.5, 3.5]]]

        mock_pdb._SSTrajectory__get_all_proteins.return_value = mock_proteins
        mock_pdb.proteinTrajectoryList = [mock_proteins]

        MockSSTrajectory.return_value = mock_pdb

        result = osmofold_local.get_chain_info("dummy.pdb")
        
        self.assertEqual(result["Chain 1"], ("ACDE", [1.1, 2.2, 3.3, 4.4], [0.5, 1.5, 2.5, 3.5]))
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_valid_pdb_multiple_chains(self, MockSSTrajectory):
        mock_pdb = MagicMock()
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = "ACDEFGH"
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_all_SASA.side_effect = [[list(range(1, 8))], [list(range(8, 15))]]
        mock_chain1 = MagicMock()
        mock_chain1.get_amino_acid_sequence.return_value = "ACD"
        mock_chain2 = MagicMock()
        mock_chain2.get_amino_acid_sequence.return_value = "EFGH"
        mock_pdb.proteinTrajectoryList = [mock_chain1, mock_chain2]
        MockSSTrajectory.return_value = mock_pdb
        
        result = osmofold_local.get_chain_info("dummy.pdb")
        self.assertEqual(result["Chain 1"], ("ACD", [1, 2, 3], [8, 9, 10]))
        self.assertEqual(result["Chain 2"], ("EFGH", [4, 5, 6, 7], [11, 12, 13, 14]))
        self.assertEqual(result["All"], ("ACDEFGH", list(range(1, 8)), list(range(8, 15))))
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_empty_pdb(self, MockSSTrajectory):
        mock_pdb = MagicMock()
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = ""
        mock_pdb._SSTrajectory__get_all_proteins.return_value.get_all_SASA.side_effect = [[], []]
        mock_pdb.proteinTrajectoryList = []
        MockSSTrajectory.return_value = mock_pdb
        
        with self.assertRaises(Exception) as context:
            osmofold_local.get_chain_info("invalid.pdb")
        self.assertEqual(str(context.exception), "list index out of range")
    
    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_invalid_pdb(self, MockSSTrajectory):
        MockSSTrajectory.side_effect = Exception("Invalid PDB format")
        
        with self.assertRaises(Exception) as context:
            osmofold_local.get_chain_info("invalid.pdb")
        self.assertEqual(str(context.exception), "Invalid PDB format")

class TestProteinUnfoldedDG(unittest.TestCase):
      
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    def test_single_osmolyte(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = (1.0 * 5 + 0.5 * 5)
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    def test_multiple_osmolytes(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5), ([0.8] * 5, [0.3] * 5)]
        pdb = "test.pdb"
        osmolytes = ["urea", "trehalose"]
        expected_dG = {
            "urea": (1.0 * 5 + 0.5 * 5),
            "trehalose": (0.8 * 5 + 0.3 * 5)
        }
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_local.extract_sequences_by_chains")
    @patch("osmofold.osmofold_local.get_tfe")
    def test_split_chains(self, mock_get_tfe, mock_extract_sequences_by_chains):
        mock_extract_sequences_by_chains.return_value = ["ABCDE", "FGHIJ"]
        
        # Three return values: Two for chains, one for the full sequence
        mock_get_tfe.side_effect = [
            ([1.0] * 5, [0.5] * 5),  # Chain 1
            ([0.8] * 5, [0.3] * 5),  # Chain 2
            ([0.9] * 10, [0.4] * 10)  # Full sequence (sum of chains)
        ]
        
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = {
            "Chain 1": {"urea": (1.0 * 5 + 0.5 * 5)},
            "Chain 2": {"urea": (0.8 * 5 + 0.3 * 5)},
            "All": {"urea": (0.9 * 10 + 0.4 * 10)}  # Full sequence calculation
        }

        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes, split_chains=True)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    def test_custom_tfe(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.return_value = ([0.5] * 5, [0.2] * 5)
        pdb = "test.pdb"
        osmolytes = "tfe"
        custom_tfe = {"tfe": 1.0}
        expected_dG = (0.5 * 5 + 0.2 * 5)
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes, custom_tfe=custom_tfe)
        self.assertEqual(result, {"tfe": expected_dG})
    
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    def test_no_osmolytes(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        pdb = "test.pdb"
        osmolytes = []
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes)
        self.assertEqual(result, {})  # Expecting an empty dictionary

class TestProteinFoldedDG(unittest.TestCase):

    @patch("osmofold.osmofold_local.get_pdb_info")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_single_osmolyte(self, mock_sasa_to_rasa, mock_get_tfe, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ABCDE", [1.0] * 5, [0.5] * 5)
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = (1.0 * 0.8 + 0.5 * 0.5) * 5
        
        result = osmofold_local.protein_folded_dG(pdb, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_local.get_pdb_info")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_multiple_osmolytes(self, mock_sasa_to_rasa, mock_get_tfe, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ABCDE", [1.0] * 5, [0.5] * 5)
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5), ([0.8] * 5, [0.3] * 5)]
        pdb = "test.pdb"
        osmolytes = ["urea", "tmao"]
        expected_dG = {
            "urea": (1.0 * 0.8 + 0.5 * 0.5) * 5,
            "tmao": (0.8 * 0.8 + 0.3 * 0.5) * 5
        }
        
        result = osmofold_local.protein_folded_dG(pdb, osmolytes)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_local.get_chain_info")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_split_chains(self, mock_sasa_to_rasa, mock_get_tfe, mock_get_chain_info):
        mock_get_chain_info.return_value = {
            "Chain 1": ("ABCDE", [1.0] * 5, [0.5] * 5), 
            "Chain 2": ("FGHIJ", [0.8] * 5, [0.3] * 5),
            "All": ("ABCDEFGHIJ", [1.0] * 5 + [0.8] * 5, [0.5] * 5 + [0.3] * 5)
        }
        mock_sasa_to_rasa.side_effect = [([0.8] * 5, [0.5] * 5), ([0.7] * 5, [0.4] * 5), ([1.0] * 10, [0.5] * 10)]

        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5),  # Chain 1
                            ([0.8] * 5, [0.3] * 5),  # Chain 2
                            ([1.0] * 10, [0.5] * 10)]  # "All"

        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = {
            "Chain 1": {"urea": (1.0 * 0.8 + 0.5 * 0.5) * 5},
            "Chain 2": {"urea": (0.8 * 0.7 + 0.3 * 0.4) * 5}
        }
        
        result = osmofold_local.protein_folded_dG(pdb, osmolytes, split_chains=True)
        print(result)
        
        for chain in expected_dG:
            for osmolyte in expected_dG[chain]:
                self.assertAlmostEqual(result[chain][osmolyte], expected_dG[chain][osmolyte], places=4)


    @patch("osmofold.osmofold_local.get_pdb_info")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_no_osmolytes(self, mock_sasa_to_rasa, mock_get_tfe, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ABCDE", [1.0] * 5, [0.5] * 5)
        pdb = "test.pdb"
        osmolytes = []
        
        result = osmofold_local.protein_folded_dG(pdb, osmolytes)
        self.assertEqual(result, {})  # Expecting an empty dictionary

class TestProteinDDGFolding(unittest.TestCase):

    @patch("osmofold.osmofold_local.protein_folded_dG")
    @patch("osmofold.osmofold_local.protein_unfolded_dG")
    def test_single_osmolyte(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0}
        mock_unfolded_dG.return_value = {"urea": -5.0}
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_ddG = {"urea": -5.0}
        
        result = osmofold_local.protein_ddG_folding(pdb, osmolytes)
        self.assertEqual(result, expected_ddG)
    
    @patch("osmofold.osmofold_local.protein_folded_dG")
    @patch("osmofold.osmofold_local.protein_unfolded_dG")
    def test_multiple_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0, "tmao": -12.0}
        mock_unfolded_dG.return_value = {"urea": -5.0, "tmao": -6.0}
        pdb = "test.pdb"
        osmolytes = ["urea", "tmao"]
        expected_ddG = {"urea": -5.0, "tmao": -6.0}
        
        result = osmofold_local.protein_ddG_folding(pdb, osmolytes)
        self.assertEqual(result, expected_ddG)
    
    @patch("osmofold.osmofold_local.protein_folded_dG")
    @patch("osmofold.osmofold_local.protein_unfolded_dG")
    def test_triplet_output(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0}
        mock_unfolded_dG.return_value = {"urea": -5.0}
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_result = {"urea": (-10.0, -5.0, -5.0)}
        
        result = osmofold_local.protein_ddG_folding(pdb, osmolytes, triplet=True)
        self.assertEqual(result, expected_result)
    
    @patch("osmofold.osmofold_local.protein_folded_dG")
    @patch("osmofold.osmofold_local.protein_unfolded_dG")
    def test_split_chains(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"Chain 1": {"urea": -10.0}, "Chain 2": {"urea": -12.0}}
        mock_unfolded_dG.return_value = {"Chain 1": {"urea": -5.0}, "Chain 2": {"urea": -6.0}}
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_ddG = {
            "Chain 1": {"urea": -5.0},
            "Chain 2": {"urea": -6.0}
        }
        
        result = osmofold_local.protein_ddG_folding(pdb, osmolytes, split_chains=True)
        self.assertEqual(result, expected_ddG)
    
    @patch("osmofold.osmofold_local.protein_folded_dG")
    @patch("osmofold.osmofold_local.protein_unfolded_dG")
    def test_no_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        pdb = "test.pdb"
        osmolytes = []
        
        result = osmofold_local.protein_ddG_folding(pdb, osmolytes)
        self.assertEqual(result, {})  # Expecting an empty dictionary

if __name__ == "__main__":
    unittest.main()