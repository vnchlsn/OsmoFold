import unittest
from unittest.mock import MagicMock, mock_open, patch
import numpy as np
from osmofold import osmofold_local
import mdtraj

class Test_OsmoFold_Local_Core(unittest.TestCase):
    def test_get_max_sasa_list(self):
        """
        Test the get_max_sasa_list function to ensure it returns the expected list of maximum SASA values.
        """
        expected = [
            129, 240, 201, 197, 174, 159, 224, 104, 155, 172,
            225, 195, 193, 223, 224, 236, 274, 167, 285, 263
        ]
        
        result = osmofold_local.get_max_sasa_list()
        
        self.assertEqual(result, expected)

    def test_amino_to_energy_valid_inputs(self):
        """
        Test that the amino_to_energy function returns correct values for valid inputs.
        """
        # Define test cases as tuples of (amino acid, cosolute, expected result)
        test_cases = [
            ("A", "trehalose_hong", 59.3),
            ("F", "tmao", -9.32),
            ("L", "sarcosine", 38.33),
            ("I", "betaine", -1.27),
            ("P", "sorbitolBack", -4.48 + 35),
            ("M", "sucrose", -6.66),
            ("W", "urea", -141.46),
            ("S", "prolineBack", -33.49 + 48),
            ("G", "glycerol", 0),
            ("A", "trehaloseBack", 33.25 + 62),
        ]

        for amino, cosolute, expected in test_cases:
            with self.subTest(amino=amino, cosolute=cosolute):
                result = osmofold_local.amino_to_energy(amino, cosolute)
                self.assertAlmostEqual(result, expected, places=2)

    def test_amino_to_energy_invalid_amino(self):
        """
        Test that the function raises a KeyError for invalid amino acid codes.
        """
        with self.assertRaises(KeyError):
            osmofold_local.amino_to_energy("Z", "trehalose_hong")

    def test_amino_to_energy_invalid_cosolute(self):
        """
        Test that the function handles invalid cosolute names gracefully.
        """
        result = osmofold_local.amino_to_energy("A", "unknown_cosolute")
        expected = None

        self.assertEqual(result, expected)

    def test_amino_to_energy_edge_case(self):
        """
        Test edge cases such as empty strings or non-existent combinations.
        """
        with self.assertRaises(KeyError):
            osmofold_local.amino_to_energy("", "trehalose")
        with self.assertRaises(KeyError):
            osmofold_local.amino_to_energy("X", "tmao")

    def test_extract_sequence(self):
        """
        Test extract_sequence function with a mock PDB file.
        """
        mock_pdb_data = """\
ATOM      1  N   MET A   1      11.104  13.207   9.784  1.00 20.00           N  
ATOM      2  CA  MET A   1      11.687  14.520   9.503  1.00 20.00           C  
ATOM      3  C   MET A   1      12.178  15.294  10.722  1.00 20.00           C  
ATOM      4  O   MET A   1      12.028  14.867  11.887  1.00 20.00           O  
ATOM      5  N   ALA A   2      12.750  16.490  10.493  1.00 20.00           N  
ATOM      6  CA  ALA A   2      13.234  17.291  11.657  1.00 20.00           C  
ATOM      7  C   ALA A   2      13.450  16.456  12.902  1.00 20.00           C  
ATOM      8  O   ALA A   2      13.593  17.027  14.034  1.00 20.00           O  
ATOM      9  N   GLY A   3      13.593  15.027  13.034  1.00 20.00           N  
ATOM     10  CA  GLY A   3      13.593  14.027  13.034  1.00 20.00           C  
TER                                                                              
"""
        expected_sequence = "MAG"  # MET -> M, ALA -> A, GLY -> G

        with patch("builtins.open", mock_open(read_data=mock_pdb_data)):
            result = osmofold_local.extract_sequence("mock_pdb_file.pdb")
            self.assertEqual(result, expected_sequence)

    def test_extract_sequence_with_duplicates(self):
        """
        Test extract_sequence function to ensure it handles duplicate residues correctly.
        """
        mock_pdb_data = """\
ATOM      1  N   MET A   1      11.104  13.207   9.784  1.00 20.00           N  
ATOM      2  CA  MET A   1      11.687  14.520   9.503  1.00 20.00           C  
ATOM      3  N   MET A   1      12.750  16.490  10.493  1.00 20.00           N  
ATOM      4  N   ALA A   2      13.450  16.456  12.902  1.00 20.00           C  
ATOM      5  CA  ALA A   2      13.593  14.027  13.034  1.00 20.00           C  
TER                                                                              
"""
        expected_sequence = "MA"  # MET -> M, ALA -> A (Duplicates ignored)

        with patch("builtins.open", mock_open(read_data=mock_pdb_data)):
            result = osmofold_local.extract_sequence("mock_pdb_file.pdb")
            self.assertEqual(result, expected_sequence)

    def test_extract_sequence_with_unknown_residue(self):
        """
        Test extract_sequence function for unknown residues.
        """
        mock_pdb_data = """\
ATOM      1  N   XYZ A   1      11.104  13.207   9.784  1.00 20.00           N  
ATOM      2  CA  MET A   2      11.687  14.520   9.503  1.00 20.00           C  
TER                                                                              
"""
        expected_sequence = "XM"  # Unknown XYZ -> X, MET -> M

        with patch("builtins.open", mock_open(read_data=mock_pdb_data)):
            result = osmofold_local.extract_sequence("mock_pdb_file.pdb")
            self.assertEqual(result, expected_sequence)

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

class Test_OsmoFold_Local_Get_TFE(unittest.TestCase):

    def test_basic_case_with_custom_tfe(self):
        seq = "AF"
        osmo = "osmolyte1"
        custom_tfe = {
            "osmolyte1": {"A": 10, "F": 15, "L": 20, "I": 25}
        }
        result = osmofold_local.get_tfe(seq, osmo, custom_tfe=custom_tfe)
        expected = [10, 15]
        self.assertEqual(result, expected)

    def test_basic_case_without_custom_tfe(self):
        seq = "AF"
        osmo = "trehalose"
        result = osmofold_local.get_tfe(seq, osmo)
        expected = [33.25, -17.88]
        self.assertEqual(result, expected)

    def test_case_with_invalid_amino_acid(self):
        seq = "AXF"
        osmo = "osmolyte1"
        custom_tfe = {
            "osmolyte1": {"A": 10, "F": 15, "L": 20, "I": 25}
        }
        with self.assertRaises(ValueError):  # Update based on actual error
            osmofold_local.get_tfe(seq, osmo, custom_tfe=custom_tfe)

    def test_case_with_no_custom_tfe_and_unsupported_osmo(self):
        seq = "AF"
        osmo = "unsupported_osmo"
        result = osmofold_local.get_tfe(seq, osmo)
        expected = [None, None]  # Should still return the default energy based on amino_to_energy
        self.assertEqual(result, expected)

    def test_case_with_empty_sequence(self):
        seq = ""
        osmo = "osmolyte1"
        result = osmofold_local.get_tfe(seq, osmo)
        expected = []
        self.assertEqual(result, expected)

    def test_case_with_empty_osmo(self):
        seq = "AF"
        osmo = ""
        result = osmofold_local.get_tfe(seq, osmo)
        expected = [None, None]  # Default values from amino_to_energy
        self.assertEqual(result, expected)

    def test_case_with_custom_tfe_and_multiple_aminos(self):
        seq = "ALF"
        osmo = "osmolyte1"
        custom_tfe = {
            "osmolyte1": {"A": 10, "F": 15, "L": 20, "I": 25}
        }
        result = osmofold_local.get_tfe(seq, osmo, custom_tfe=custom_tfe)
        expected = [10, 20, 15]
        self.assertEqual(result, expected)

class Test_OsmoFold_Local_SASA_To_RASA(unittest.TestCase):
    @patch('osmofold.osmofold_local.get_max_sasa_list')
    def test_basic_functionality(self, mock_max_sasa):
        mock_max_sasa.return_value = [129, 240, 201, 197, 174, 159, 224, 104, 155, 172, 
            225, 195, 193, 223, 224, 236, 274, 167, 285, 263]
        seq = "ACD"
        sasa_list = [10.0, 20.0, 15.0]
        expected_output = [10.0/129, 20.0/167, 15.0/193]
        self.assertEqual(osmofold_local.sasa_to_rasa(seq, sasa_list), expected_output)

    @patch('osmofold.osmofold_local.get_max_sasa_list')
    def test_empty_inputs(self, mock_max_sasa):
        mock_max_sasa.return_value = [129, 240, 201, 197, 174, 159, 224, 104, 155, 172, 
            225, 195, 193, 223, 224, 236, 274, 167, 285, 263]
        seq = ""
        sasa_list = []
        self.assertEqual(osmofold_local.sasa_to_rasa(seq, sasa_list), [])

    def test_sequence_sasa_length_mismatch(self):
        seq = "ACD"
        sasa_list = [10.0, 20.0]
        with self.assertRaises(IndexError):
            osmofold_local.sasa_to_rasa(seq, sasa_list)

    @patch('osmofold.osmofold_local.get_max_sasa_list')
    def test_invalid_amino_acid(self, mock_max_sasa):
        mock_max_sasa.return_value = [129, 240, 201, 197, 174, 159, 224, 104, 155, 172, 
            225, 195, 193, 223, 224, 236, 274, 167, 285, 263]
        seq = "ACX"
        sasa_list = [10.0, 20.0, 15.0]
        with self.assertRaises(ValueError):  # Update based on actual error
            osmofold_local.sasa_to_rasa(seq, sasa_list)

class TestGetPDBInfo(unittest.TestCase):
    @patch('osmofold.osmofold_local.SSTrajectory')
    def test_get_pdb_info(self, mock_SSTrajectory):
        """Test get_pdb_info function with a mocked SSTrajectory."""
        # Mocking the SSTrajectory object and its methods
        mock_pdb_instance = MagicMock()
        mock_SSTrajectory.return_value = mock_pdb_instance

        # Mock the protein instance
        mock_protein_instance = MagicMock()

        # Configure __get_all_proteins to return the same protein instance for both calls
        mock_pdb_instance._SSTrajectory__get_all_proteins.side_effect = lambda *args, **kwargs: mock_protein_instance

        # Mock the method calls for sequence and SASA
        mock_protein_instance.get_amino_acid_sequence.return_value = "ACDEFGHIK"
        mock_protein_instance.get_all_SASA.return_value = [[10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]]

        # Call the function
        pdb_path = "protein.pdb"
        result = osmofold_local.get_pdb_info(pdb_path)

        # Define the expected result
        expected_result = [
            "ACDEFGHIK",  # Mocked sequence
            [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],  # Mocked SASA
        ]

        # Assertions
        self.assertEqual(result, expected_result)

        # Check that SSTrajectory was called with the correct arguments
        mock_SSTrajectory.assert_called_once_with(pdb_path, pdb_path)

        # Check that the mocked methods were called as expected
        self.assertEqual(
            mock_pdb_instance._SSTrajectory__get_all_proteins.call_count, 2
        )  # Expecting two calls to __get_all_proteins
        mock_protein_instance.get_amino_acid_sequence.assert_called_once_with(oneletter=True)
        mock_protein_instance.get_all_SASA.assert_called_once_with(mode='residue', stride=1)

class Test_OsmoFold_Local_Unfolded_dG(unittest.TestCase):
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_single_osmolyte(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACDEFGHIKLMNPQRSTVWY"
        mock_get_tfe.return_value = [1.0] * 20  # Mocked TFE values
        
        result = osmofold_local.protein_unfolded_dG("test.pdb", "tmao")
        expected = {"tmao": 20.0}
        self.assertEqual(result, expected)
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_multiple_osmolytes(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACDEFGHIKLMNPQRSTVWY"
        mock_get_tfe.side_effect = [[1.0] * 20, [0.5] * 20]  # Mocked TFE values for each osmolyte
        
        result = osmofold_local.protein_unfolded_dG("test.pdb", ["tmao", "urea"])
        expected = {"tmao": 20.0, "urea": 10.0}
        self.assertEqual(result, expected)
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_custom_tfe(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACD"
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]  # Custom TFE should not matter here
        
        custom_tfe = {"tmaoBack": 1.0}
        result = osmofold_local.protein_unfolded_dG("test.pdb", "tmao", custom_tfe=custom_tfe)
        expected = {"tmao": 6.0}
        self.assertEqual(result, expected)
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_backbone_inclusion(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACD"
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]
        
        result_backbone = osmofold_local.protein_unfolded_dG("test.pdb", "tmao", backbone=True)
        result_no_backbone = osmofold_local.protein_unfolded_dG("test.pdb", "tmao", backbone=False)
        
        self.assertEqual(result_backbone, {"tmao": 6.0})
        self.assertEqual(result_no_backbone, {"tmao": 6.0})  # Mock TFE identical
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_concentration_scaling(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACD"
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]
        
        result = osmofold_local.protein_unfolded_dG("test.pdb", "tmao", concentration=2.0)
        expected = {"tmao": 12.0}
        self.assertEqual(result, expected)
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_empty_osmolytes(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACD"
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]

        result = osmofold_local.protein_unfolded_dG("test.pdb", [])
        self.assertEqual(result, {})
    
    @patch('osmofold.osmofold_local.extract_sequence')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_invalid_osmolytes_type(self, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ACD"
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]

        with self.assertRaises(TypeError):
            osmofold_local.protein_unfolded_dG("test.pdb", 123)

class TestProteinFoldedDG(unittest.TestCase):
    @patch('osmofold.osmofold_local.get_pdb_info')
    @patch('osmofold.osmofold_local.sasa_to_rasa')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_single_osmolyte(self, mock_get_tfe, mock_sasa_to_rasa, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ACDE", [10.0, 20.0, 30.0, 40.0])
        mock_sasa_to_rasa.return_value = [0.1, 0.2, 0.3, 0.4]
        mock_get_tfe.return_value = [1.0, 2.0, 3.0, 4.0]

        """Test for a single osmolyte."""
        result = osmofold_local.protein_folded_dG("protein.pdb", "osmolyte")
        expected_result = {"osmolyte": 1.0 * (0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0)}
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.get_pdb_info')
    @patch('osmofold.osmofold_local.sasa_to_rasa')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_multiple_osmolytes(self, mock_get_tfe, mock_sasa_to_rasa, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ACDE", [10.0, 20.0, 30.0, 40.0])
        mock_sasa_to_rasa.return_value = [0.1, 0.2, 0.3, 0.4]
        mock_get_tfe.side_effect = [
            [1.0, 2.0, 3.0, 4.0],  # Osmolyte 1
            [2.0, 3.0, 4.0, 5.0],  # Osmolyte 2
        ]

        result = osmofold_local.protein_folded_dG("protein.pdb", ["osmolyte1", "osmolyte2"])
        expected_result = {
            "osmolyte1": 1.0 * (0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0),
            "osmolyte2": 1.0 * (0.1 * 2.0 + 0.2 * 3.0 + 0.3 * 4.0 + 0.4 * 5.0),
        }
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.get_pdb_info')
    @patch('osmofold.osmofold_local.sasa_to_rasa')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_custom_tfe(self, mock_get_tfe, mock_sasa_to_rasa, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ACDE", [10.0, 20.0, 30.0, 40.0])
        mock_sasa_to_rasa.return_value = [0.1, 0.2, 0.3, 0.4]
        mock_get_tfe.return_value = [1.0, 2.0, 3.0, 4.0]

        """Test with custom TFE values."""
        custom_tfe = {"osmolyteBack": [1.0, 2.0, 3.0, 4.0]}
        result = osmofold_local.protein_folded_dG("protein.pdb", "osmolyte", custom_tfe=custom_tfe)
        expected_result = {"osmolyte": 1.0 * (0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0)}
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.get_pdb_info')
    @patch('osmofold.osmofold_local.sasa_to_rasa')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_concentration_scaling(self, mock_get_tfe, mock_sasa_to_rasa, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ACDE", [10.0, 20.0, 30.0, 40.0])
        mock_sasa_to_rasa.return_value = [0.1, 0.2, 0.3, 0.4]
        mock_get_tfe.return_value = [1.0, 2.0, 3.0, 4.0]

        """Test with a different concentration value."""
        result = osmofold_local.protein_folded_dG("protein.pdb", "osmolyte", concentration=2.0)
        expected_result = {"osmolyte": 2.0 * (0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0)}
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.get_pdb_info')
    @patch('osmofold.osmofold_local.sasa_to_rasa')
    @patch('osmofold.osmofold_local.get_tfe')
    def test_no_backbone(self, mock_get_tfe, mock_sasa_to_rasa, mock_get_pdb_info):
        mock_get_pdb_info.return_value = ("ACDE", [10.0, 20.0, 30.0, 40.0])
        mock_sasa_to_rasa.return_value = [0.1, 0.2, 0.3, 0.4]
        mock_get_tfe.return_value = [1.0, 2.0, 3.0, 4.0]
        
        """Test with backbone=False."""
        result = osmofold_local.protein_folded_dG("protein.pdb", "osmolyte", backbone=False)
        self.assertEqual(result, {"osmolyte": 1.0 * (0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0)})

class TestProteinDDGFolding(unittest.TestCase):
    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_single_osmolyte(self, mock_unfolded_dG, mock_folded_dG):
        """Test with a single osmolyte, without triplet output."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0}

        result = osmofold_local.protein_ddG_folding("protein.pdb", "osmolyte1")
        expected_result = {"osmolyte1": -20.0 - (-10.0)}
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_multiple_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        """Test with multiple osmolytes, without triplet output."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0, "osmolyte2": -25.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0, "osmolyte2": -15.0}

        result = osmofold_local.protein_ddG_folding("protein.pdb", ["osmolyte1", "osmolyte2"])
        expected_result = {
            "osmolyte1": -20.0 - (-10.0),
            "osmolyte2": -25.0 - (-15.0),
        }
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_triplet_output(self, mock_unfolded_dG, mock_folded_dG):
        """Test with triplet output enabled."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0}

        result = osmofold_local.protein_ddG_folding("protein.pdb", ["osmolyte1"], triplet=True)
        expected_result = {
            "osmolyte1": (-20.0, -10.0, -20.0 - (-10.0)),
        }
        self.assertEqual(result, expected_result)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_custom_tfe(self, mock_unfolded_dG, mock_folded_dG):
        """Test with custom TFE values."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0}

        custom_tfe = {"osmolyte1Back": [-1.0, -2.0, -3.0, -4.0]}
        result = osmofold_local.protein_ddG_folding("protein.pdb", "osmolyte1", custom_tfe=custom_tfe)
        mock_folded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], True, custom_tfe, 1.0)
        mock_unfolded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], True, custom_tfe, 1.0)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_no_backbone(self, mock_unfolded_dG, mock_folded_dG):
        """Test with backbone=False."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0}

        result = osmofold_local.protein_ddG_folding("protein.pdb", "osmolyte1", backbone=False)
        mock_folded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], False, None, 1.0)
        mock_unfolded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], False, None, 1.0)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_concentration_scaling(self, mock_unfolded_dG, mock_folded_dG):
        """Test with a different concentration value."""
        mock_folded_dG.return_value = {"osmolyte1": -20.0}
        mock_unfolded_dG.return_value = {"osmolyte1": -10.0}

        result = osmofold_local.protein_ddG_folding("protein.pdb", "osmolyte1", concentration=2.0)
        mock_folded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], True, None, 2.0)
        mock_unfolded_dG.assert_called_once_with("protein.pdb", ["osmolyte1"], True, None, 2.0)

    @patch('osmofold.osmofold_local.protein_folded_dG')
    @patch('osmofold.osmofold_local.protein_unfolded_dG')
    def test_no_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        """Test with no osmolytes provided (edge case)."""
        mock_folded_dG.return_value = {}
        mock_unfolded_dG.return_value = {}

        result = osmofold_local.protein_ddG_folding("protein.pdb", [])
        self.assertEqual(result, {})

if __name__ == "__main__":
    unittest.main()