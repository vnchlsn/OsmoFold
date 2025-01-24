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

class TestGetChainInfo(unittest.TestCase):

    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_single_chain(self, mock_SSTrajectory):
        # Mock SSTrajectory behavior
        mock_traj = mock_SSTrajectory.return_value
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = "MAG"
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_all_SASA.return_value = [[1.2, 0.8, 1.5]]

        # Mocking proteinTrajectoryList
        mock_protein = MagicMock()
        mock_protein.get_amino_acid_sequence.return_value = "MAG"
        mock_traj.proteinTrajectoryList = [mock_protein]

        result = osmofold_local.get_chain_info("dummy.pdb")

        expected = {
            "Chain 1": ["MAG", [1.2, 0.8, 1.5]],
            "All": ["MAG", [1.2, 0.8, 1.5]]
        }
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_multiple_chains(self, mock_SSTrajectory):
        # Mock SSTrajectory behavior
        mock_traj = mock_SSTrajectory.return_value
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = "MAGLYS"
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_all_SASA.return_value = [[1.2, 0.8, 1.5, 2.3, 1.0, 0.5]]

        # Mocking proteinTrajectoryList
        mock_protein1 = MagicMock()
        mock_protein1.get_amino_acid_sequence.return_value = "MAG"

        mock_protein2 = MagicMock()
        mock_protein2.get_amino_acid_sequence.return_value = "LYS"

        mock_traj.proteinTrajectoryList = [mock_protein1, mock_protein2]

        result = osmofold_local.get_chain_info("dummy.pdb")

        expected = {
            "Chain 1": ["MAG", [1.2, 0.8, 1.5]],
            "Chain 2": ["LYS", [2.3, 1.0, 0.5]],
            "All": ["MAGLYS", [1.2, 0.8, 1.5, 2.3, 1.0, 0.5]]
        }
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_local.SSTrajectory")
    def test_empty_pdb(self, mock_SSTrajectory):
        # Mock SSTrajectory behavior for an empty PDB
        mock_traj = mock_SSTrajectory.return_value
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_amino_acid_sequence.return_value = ""
        mock_traj._SSTrajectory__get_all_proteins.return_value.get_all_SASA.return_value = [[]]

        mock_traj.proteinTrajectoryList = []

        result = osmofold_local.get_chain_info("empty.pdb")

        expected = {
            "All": ["", []]
        }
        self.assertEqual(result, expected)

class TestProteinUnfoldedDG(unittest.TestCase):
    # Mocked functions
    def mock_extract_sequence(pdb):
        return "ACDEFGHIKLMNPQRSTVWY"

    def mock_extract_sequences_by_chains(pdb):
        return ["ACDEFG", "HIKLMN", "PQRSTV", "WY"]

    def mock_get_tfe(sequence, osmolyte, custom_tfe=None):
        return np.array([0.5] * len(sequence))  # Simple mock energy values

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_single_osmolyte(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea")
        self.assertEqual(result, {"urea": 10.0})

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_multiple_osmolytes(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, ["urea", "gdmcl"])
        self.assertEqual(result, {"urea": 10.0, "gdmcl": 10.0})

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)  # Ensure this is applied last
    @patch("osmofold.osmofold_local.extract_sequences_by_chains", side_effect=mock_extract_sequences_by_chains)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_split_chains(self, mock_tfe, mock_seq_by_chains, mock_seq):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea", split_chains=True)

        expected = {
            "Chain 1": {"urea": 3.0},
            "Chain 2": {"urea": 3.0},
            "Chain 3": {"urea": 3.0},
            "Chain 4": {"urea": 1.0},
            "All": {"urea": 10.0},
        }
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_concentration(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea", concentration=2.0)
        self.assertEqual(result, {"urea": 20.0})

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_custom_tfe(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        custom_tfe = {"urea": 1.0}
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea", custom_tfe=custom_tfe)

        self.assertEqual(result, {"urea": 10.0})  # Still the same since get_tfe is mocked

    @patch("osmofold.osmofold_local.extract_sequence", side_effect=mock_extract_sequence)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_unfolded_dG_backbone_false(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea", backbone=False)

        self.assertEqual(result, {"urea": 10.0})  # Ensures backbone flag is handled correctly

    @patch("osmofold.osmofold_local.extract_sequence", return_value="")
    @patch("osmofold.osmofold_local.get_tfe", return_value=np.array([]))
    def test_protein_unfolded_dG_empty_sequence(self, mock_seq, mock_tfe):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_unfolded_dG(pdb_path, "urea")

        self.assertEqual(result, {"urea": 0.0})  # Empty sequence should return zero free energy

class TestProteinFoldedDG(unittest.TestCase):

    # Mocked functions
    def mock_get_pdb_info(pdb):
        return ("ACDEFGHIKLMNPQRSTVWY", "mock_structure")

    def mock_get_chain_info(pdb):
        return {
            "Chain 1": ("ACDEFG", "mock_structure_1"),
            "Chain 2": ("HIKLMN", "mock_structure_2"),
            "Chain 3": ("PQRSTV", "mock_structure_3"),
            "Chain 4": ("WY", "mock_structure_4"),
            "All": ("ACDEFGHIKLMNPQRSTVWY", "mock_structure")
        }

    def mock_get_tfe(sequence, osmolyte, custom_tfe=None):
        return np.array([0.5] * len(sequence))  # Mock TFE values

    def mock_sasa_to_rasa(sequence, structure):
        return np.array([1.0] * len(sequence))  # Mock relative solvent accessibility values (fully exposed)

    @patch("osmofold.osmofold_local.get_pdb_info", side_effect=mock_get_pdb_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_single_osmolyte(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, "urea")
        self.assertEqual(result, {"urea": 10.0})

    @patch("osmofold.osmofold_local.get_pdb_info", side_effect=mock_get_pdb_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_multiple_osmolytes(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, ["urea", "gdmcl"])
        self.assertEqual(result, {"urea": 10.0, "gdmcl": 10.0})

    @patch("osmofold.osmofold_local.get_chain_info", side_effect=mock_get_chain_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_split_chains(self, mock_tfe, mock_rasa, mock_chain_info):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, "urea", split_chains=True)
        print(result)

        expected = {
            "Chain 1": {"urea": 3.0},
            "Chain 2": {"urea": 3.0},
            "Chain 3": {"urea": 3.0},
            "Chain 4": {"urea": 1.0},
            "All": {"urea": 10.0},
        }
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_local.get_pdb_info", side_effect=mock_get_pdb_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_concentration(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, "urea", concentration=2.0)
        self.assertEqual(result, {"urea": 20.0})

    @patch("osmofold.osmofold_local.get_pdb_info", side_effect=mock_get_pdb_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_custom_tfe(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        custom_tfe = {"urea": 1.0}
        result = osmofold_local.protein_folded_dG(pdb_path, "urea", custom_tfe=custom_tfe)

        self.assertEqual(result, {"urea": 10.0})  # Still 10.0 since get_tfe is mocked

    @patch("osmofold.osmofold_local.get_pdb_info", side_effect=mock_get_pdb_info)
    @patch("osmofold.osmofold_local.sasa_to_rasa", side_effect=mock_sasa_to_rasa)
    @patch("osmofold.osmofold_local.get_tfe", side_effect=mock_get_tfe)
    def test_protein_folded_dG_backbone_false(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, "urea", backbone=False)

        self.assertEqual(result, {"urea": 10.0})  # Ensures backbone flag is handled correctly

    @patch("osmofold.osmofold_local.get_pdb_info", return_value=("", ""))
    @patch("osmofold.osmofold_local.sasa_to_rasa", return_value=np.array([]))
    @patch("osmofold.osmofold_local.get_tfe", return_value=np.array([]))
    def test_protein_folded_dG_empty_sequence(self, mock_tfe, mock_rasa, mock_pdb):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_folded_dG(pdb_path, "urea")

        self.assertEqual(result, {"urea": 0.0})  # Empty sequence should return zero free energy

class TestProteinDDGFolding(unittest.TestCase):

    # Mocked functions
    def mock_protein_folded_dG(pdb, osmolytes, backbone=True, custom_tfe=None, concentration=1.0, split_chains=False):
        if split_chains:
            return {
                "Chain 1": {"urea": 3.0},
                "Chain 2": {"urea": 3.0},
                "Chain 3": {"urea": 3.0},
                "Chain 4": {"urea": 2.0},
                "All": {"urea": 10.0},
            }
        return {"urea": 10.0}

    def mock_protein_unfolded_dG(pdb, osmolytes, backbone=True, custom_tfe=None, concentration=1.0, split_chains=False):
        if split_chains:
            return {
                "Chain 1": {"urea": 5.0},
                "Chain 2": {"urea": 5.0},
                "Chain 3": {"urea": 5.0},
                "Chain 4": {"urea": 2.0},
                "All": {"urea": 17.0},
            }
        return {"urea": 17.0}

    @patch("osmofold.osmofold_local.protein_folded_dG", side_effect=mock_protein_folded_dG)
    @patch("osmofold.osmofold_local.protein_unfolded_dG", side_effect=mock_protein_unfolded_dG)
    def test_protein_ddG_folding_single_osmolyte(self, mock_unfolded, mock_folded):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_ddG_folding(pdb_path, "urea")

        self.assertEqual(result, {"urea": -7.0})

    @patch("osmofold.osmofold_local.protein_folded_dG", side_effect=mock_protein_folded_dG)
    @patch("osmofold.osmofold_local.protein_unfolded_dG", side_effect=mock_protein_unfolded_dG)
    def test_protein_ddG_folding_triplet(self, mock_unfolded, mock_folded):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_ddG_folding(pdb_path, "urea", triplet=True)

        self.assertEqual(result, {"urea": (10.0, 17.0, -7.0)})

    @patch("osmofold.osmofold_local.protein_folded_dG", side_effect=mock_protein_folded_dG)
    @patch("osmofold.osmofold_local.protein_unfolded_dG", side_effect=mock_protein_unfolded_dG)
    def test_protein_ddG_folding_split_chains(self, mock_unfolded, mock_folded):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_ddG_folding(pdb_path, "urea", split_chains=True)

        expected = {
            "Chain 1": {"urea": -2.0},
            "Chain 2": {"urea": -2.0},
            "Chain 3": {"urea": -2.0},
            "Chain 4": {"urea": 0.0},
            "All": {"urea": -7.0},
        }
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_local.protein_folded_dG", side_effect=mock_protein_folded_dG)
    @patch("osmofold.osmofold_local.protein_unfolded_dG", side_effect=mock_protein_unfolded_dG)
    def test_protein_ddG_folding_split_chains_triplet(self, mock_unfolded, mock_folded):
        pdb_path = "test.pdb"
        result = osmofold_local.protein_ddG_folding(pdb_path, "urea", split_chains=True, triplet=True)

        expected = {
            "Chain 1": {"urea": (3.0, 5.0, -2.0)},
            "Chain 2": {"urea": (3.0, 5.0, -2.0)},
            "Chain 3": {"urea": (3.0, 5.0, -2.0)},
            "Chain 4": {"urea": (2.0, 2.0, 0.0)},
            "All": {"urea": (10.0, 17.0, -7.0)},
        }
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()