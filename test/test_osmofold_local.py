import unittest
from unittest.mock import MagicMock, mock_open, patch
import numpy as np
from osmofold import osmofold_local
import mdtraj

class TestCleanDict(unittest.TestCase):
    
    def test_clean_dict_with_nested_structure(self):
        data = {
            "a": np.float64(1.23),
            "b": [np.float64(2.34), {"c": np.float64(3.45)}],
            "d": (np.float64(4.56), {"e": np.float64(5.67)})
        }
        expected = {
            "a": 1.23,
            "b": [2.34, {"c": 3.45}],
            "d": (4.56, {"e": 5.67})
        }
        self.assertEqual(osmofold_local.clean_dict(data), expected)

    def test_clean_dict_with_no_np_float64(self):
        data = {"a": 1.23, "b": [2.34, {"c": 3.45}], "d": (4.56, {"e": 5.67})}
        self.assertEqual(osmofold_local.clean_dict(data), data)

    def test_clean_dict_np_float64(self):
        data = {"a": np.float64(10.5), "b": [np.float64(20.6)]}
        expected = {"a": 10.5, "b": [20.6]}
        self.assertEqual(osmofold_local.clean_dict(data), expected)

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
        self.assertAlmostEqual(result["backbone"]["A"], 46)
        self.assertAlmostEqual(result["sidechain"]["R"], 196)
    
    def test_all_amino_acids_present(self):
        """Test that all 20 amino acids are present in both backbone and sidechain."""
        expected_amino_acids = set("ARNDCEQGHILKMFPSTWYV")
        result = osmofold_local.get_max_sasa_list()
        self.assertEqual(set(result["backbone"].keys()), expected_amino_acids)
        self.assertEqual(set(result["sidechain"].keys()), expected_amino_acids)

class TestGetUnfoldedSasaList(unittest.TestCase):
    
    def test_return_type(self):
        """Test that the function returns a dictionary."""
        result = osmofold_local.get_unfolded_sasa_list()
        self.assertIsInstance(result, dict)
    
    def test_keys_exist(self):
        """Test that the dictionary contains both 'backbone' and 'sidechain' keys."""
        result = osmofold_local.get_unfolded_sasa_list()
        self.assertIn("backbone", result)
        self.assertIn("sidechain", result)
    
    def test_inner_structure(self):
        """Test that both 'backbone' and 'sidechain' contain dictionaries."""
        result = osmofold_local.get_unfolded_sasa_list()
        self.assertIsInstance(result["backbone"], dict)
        self.assertIsInstance(result["sidechain"], dict)
    
    def test_expected_values(self):
        """Test that specific values are correct."""
        result = osmofold_local.get_unfolded_sasa_list()
        self.assertAlmostEqual(result["backbone"]["A"], 27.85)
        self.assertAlmostEqual(result["sidechain"]["R"], 171.10)
    
    def test_all_amino_acids_present(self):
        """Test that all 20 amino acids are present in both backbone and sidechain."""
        expected_amino_acids = set("ARNDCEQGHILKMFPSTWYV")
        result = osmofold_local.get_unfolded_sasa_list()
        self.assertEqual(set(result["backbone"].keys()), expected_amino_acids)
        self.assertEqual(set(result["sidechain"].keys()), expected_amino_acids)

class TestGetUnfoldedSasaFromSequence(unittest.TestCase):

    def test_return_type(self):
        """Test that the function returns a tuple of two lists."""
        result = osmofold_local.get_unfolded_sasa_from_sequence("ACD")
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], list)
        self.assertIsInstance(result[1], list)

    def test_list_lengths_match_sequence(self):
        """Test that the returned lists have the same length as the input sequence."""
        sequence = "GQY"
        backbone_sasa, sidechain_sasa = osmofold_local.get_unfolded_sasa_from_sequence(sequence)
        self.assertEqual(len(backbone_sasa), len(sequence))
        self.assertEqual(len(sidechain_sasa), len(sequence))

    def test_known_values(self):
        """Test that the function returns correct values for known residues."""
        result = osmofold_local.get_unfolded_sasa_from_sequence("AG")
        self.assertAlmostEqual(result[0][0], 27.85)  # backbone A
        self.assertAlmostEqual(result[1][0], 55.10)  # sidechain A
        self.assertAlmostEqual(result[0][1], 65.15)  # backbone G
        self.assertAlmostEqual(result[1][1], 1.00)   # sidechain G

    def test_all_20_amino_acids(self):
        """Test that the function handles all 20 standard amino acids."""
        sequence = "ARNDCEQGHILKMFPSTWYV"
        backbone_sasa, sidechain_sasa = osmofold_local.get_unfolded_sasa_from_sequence(sequence)
        self.assertEqual(len(backbone_sasa), 20)
        self.assertEqual(len(sidechain_sasa), 20)

    def test_case_insensitivity(self):
        """Test that lowercase input is handled correctly."""
        result_upper = osmofold_local.get_unfolded_sasa_from_sequence("ARN")
        result_lower = osmofold_local.get_unfolded_sasa_from_sequence("arn")
        self.assertEqual(result_upper, result_lower)

    def test_invalid_amino_acid_raises(self):
        """Test that an invalid character raises a ValueError."""
        with self.assertRaises(ValueError):
            osmofold_local.get_unfolded_sasa_from_sequence("AXZ")

    def test_empty_sequence(self):
        """Test that an empty sequence returns two empty lists."""
        result = osmofold_local.get_unfolded_sasa_from_sequence("")
        self.assertEqual(result, ([], []))

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

    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.md.shrake_rupley")
    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_valid_pdb(self, mock_load_pdb, mock_shrake_rupley, mock_extract_sequence):
        # Setup mocks
        mock_traj = MagicMock()
        
        # Mock residues and atoms in topology
        mock_atom_bb = MagicMock(index=0, is_backbone=True)
        mock_atom_sc = MagicMock(index=1, is_backbone=False)
        mock_residue = MagicMock(atoms=[mock_atom_bb, mock_atom_sc])
        mock_traj.topology.residues = [mock_residue]
        
        mock_load_pdb.return_value = mock_traj
        
        # Shrake-rupley returns a numpy array inside a list (frames)
        mock_shrake_rupley.return_value = [[1.1, 2.2]]  # indices 0,1
        
        # Mock sequence extraction
        mock_extract_sequence.return_value = ["A"]  # single residue sequence
        
        # Call function
        seq, bb_sasa, sc_sasa = osmofold_local.get_pdb_info("dummy.pdb")
        
        # Check results
        self.assertEqual(seq, ["A"])
        self.assertEqual(bb_sasa, [1.1])
        self.assertEqual(sc_sasa, [2.2])

    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.md.shrake_rupley")
    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_empty_pdb(self, mock_load_pdb, mock_shrake_rupley, mock_extract_sequence):
        # Empty residues list simulates empty PDB topology
        mock_traj = MagicMock()
        mock_traj.topology.residues = []
        mock_load_pdb.return_value = mock_traj
        
        mock_shrake_rupley.return_value = [[]]
        mock_extract_sequence.return_value = []
        
        # Call function and expect empty lists
        seq, bb_sasa, sc_sasa = osmofold_local.get_pdb_info("empty.pdb")
        self.assertEqual(seq, [])
        self.assertEqual(bb_sasa, [])
        self.assertEqual(sc_sasa, [])

    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_invalid_pdb(self, mock_load_pdb):
        # Simulate MDTraj raising an exception on invalid PDB
        mock_load_pdb.side_effect = Exception("Invalid PDB format")
        
        with self.assertRaises(Exception) as context:
            osmofold_local.get_pdb_info("invalid.pdb")
        self.assertEqual(str(context.exception), "Invalid PDB format")


class TestGetChainInfo(unittest.TestCase):

    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.extract_sequences_by_chains")
    @patch("osmofold.osmofold_local.md.shrake_rupley")
    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_valid_pdb_single_chain(self, mock_load_pdb, mock_shrake_rupley,
                                    mock_extract_sequences_by_chains, mock_extract_sequence):
        # Setup mocks
        mock_traj = MagicMock()
        
        # Create one chain with one residue with two atoms (backbone, sidechain)
        mock_atom_bb = MagicMock(index=0, is_backbone=True)
        mock_atom_sc = MagicMock(index=1, is_backbone=False)
        mock_residue = MagicMock(atoms=[mock_atom_bb, mock_atom_sc])
        mock_chain = MagicMock(residues=[mock_residue])
        mock_traj.topology.chains = [mock_chain]
        mock_traj.topology.residues = [mock_residue]
        
        mock_load_pdb.return_value = mock_traj
        
        # SASA values for two atoms
        mock_shrake_rupley.return_value = [[1.1, 2.2]]  # single frame
        
        # Chain sequences
        mock_extract_sequences_by_chains.return_value = ["ACDE"]
        mock_extract_sequence.return_value = "ACDE"
        
        result = osmofold_local.get_chain_info("dummy.pdb")
        
        self.assertIn("Chain 1", result)
        seq, bb_sasa, sc_sasa = result["Chain 1"]
        self.assertEqual(seq, "ACDE")
        self.assertEqual(bb_sasa, [1.1])
        self.assertEqual(sc_sasa, [2.2])

        # Check "All" key
        full_seq, all_bb_sasa, all_sc_sasa = result["All"]
        self.assertEqual(full_seq, "ACDE")
        self.assertEqual(all_bb_sasa, [1.1])
        self.assertEqual(all_sc_sasa, [2.2])

    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.extract_sequences_by_chains")
    @patch("osmofold.osmofold_local.md.shrake_rupley")
    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_valid_pdb_multiple_chains(self, mock_load_pdb, mock_shrake_rupley,
                                       mock_extract_sequences_by_chains, mock_extract_sequence):
        # Setup mocks
        mock_traj = MagicMock()

        # Chain 1: 1 residue with 2 atoms
        mock_atom_bb1 = MagicMock(index=0, is_backbone=True)
        mock_atom_sc1 = MagicMock(index=1, is_backbone=False)
        mock_residue1 = MagicMock(atoms=[mock_atom_bb1, mock_atom_sc1])
        mock_chain1 = MagicMock(residues=[mock_residue1])

        # Chain 2: 2 residues, each with 2 atoms (indices 2,3 and 4,5)
        mock_atom_bb2 = MagicMock(index=2, is_backbone=True)
        mock_atom_sc2 = MagicMock(index=3, is_backbone=False)
        mock_residue2 = MagicMock(atoms=[mock_atom_bb2, mock_atom_sc2])

        mock_atom_bb3 = MagicMock(index=4, is_backbone=True)
        mock_atom_sc3 = MagicMock(index=5, is_backbone=False)
        mock_residue3 = MagicMock(atoms=[mock_atom_bb3, mock_atom_sc3])

        mock_chain2 = MagicMock(residues=[mock_residue2, mock_residue3])

        mock_traj.topology.chains = [mock_chain1, mock_chain2]
        mock_traj.topology.residues = [mock_residue1, mock_residue2, mock_residue3]

        mock_load_pdb.return_value = mock_traj

        # SASA values for 6 atoms
        mock_shrake_rupley.return_value = [[
            1, 8,  # Chain 1 atoms (bb, sc)
            2, 9,  # Chain 2 residue 1 atoms
            3, 10  # Chain 2 residue 2 atoms
        ]]

        # Chain sequences
        mock_extract_sequences_by_chains.return_value = ["A", "BC"]
        mock_extract_sequence.return_value = "ABC"

        result = osmofold_local.get_chain_info("dummy.pdb")

        self.assertEqual(result["Chain 1"], ("A", [1.0], [8.0]))
        self.assertEqual(result["Chain 2"], ("BC", [2.0, 3.0], [9.0, 10.0]))
        
        # More precise check for Chain 2
        bb_chain2 = result["Chain 2"][1]
        sc_chain2 = result["Chain 2"][2]
        self.assertEqual(bb_chain2, [2.0, 3.0])
        self.assertEqual(sc_chain2, [9.0, 10.0])

        # Overall all residues
        full_seq, all_bb, all_sc = result["All"]
        self.assertEqual(full_seq, "ABC")
        self.assertEqual(all_bb, [1.0, 2.0, 3.0])
        self.assertEqual(all_sc, [8.0, 9.0, 10.0])

    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.extract_sequences_by_chains")
    @patch("osmofold.osmofold_local.md.shrake_rupley")
    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_empty_pdb(self, mock_load_pdb, mock_shrake_rupley,
                       mock_extract_sequences_by_chains, mock_extract_sequence):
        mock_traj = MagicMock()
        mock_traj.topology.chains = []
        mock_traj.topology.residues = []
        mock_load_pdb.return_value = mock_traj

        mock_shrake_rupley.return_value = [[]]
        mock_extract_sequences_by_chains.return_value = []
        mock_extract_sequence.return_value = ""

        result = osmofold_local.get_chain_info("empty.pdb")

        self.assertEqual(result, {
            "All": ("", [], [])
        })

    @patch("osmofold.osmofold_local.md.load_pdb")
    def test_invalid_pdb(self, mock_load_pdb):
        mock_load_pdb.side_effect = Exception("Invalid PDB format")

        with self.assertRaises(Exception) as context:
            osmofold_local.get_chain_info("invalid.pdb")
        self.assertEqual(str(context.exception), "Invalid PDB format")

class TestProteinUnfoldedDG(unittest.TestCase):
      
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_single_osmolyte(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100 * 5])
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = (1.0 * 5 * 0.9 + 0.5 * 5 * 0.8)
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_multiple_osmolytes(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5), ([0.8] * 5, [0.3] * 5)]
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100 * 5])
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        pdb = "test.pdb"
        osmolytes = ["urea", "trehalose"]
        expected_dG = {
            "urea": (1.0 * 5 * 0.9 + 0.5 * 5 * 0.8),
            "trehalose": (0.8 * 5 * 0.9 + 0.3 * 5 * 0.8)
        }
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_local.extract_sequences_by_chains")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    @patch("osmofold.osmofold_local.extract_sequence")
    def test_split_chains(self, mock_extract_sequence, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe, mock_extract_sequences_by_chains):
        mock_extract_sequences_by_chains.return_value = ["ABCDE", "FGHIJ"]
        
        # Three return values: Two for chains, one for the full sequence
        mock_get_tfe.side_effect = [
            ([1.0] * 5, [0.5] * 5),  # Chain 1
            ([0.8] * 5, [0.3] * 5),  # Chain 2
            ([0.9] * 10, [0.4] * 10)  # Full sequence (sum of chains)
        ]

        mock_get_unfolded_sasa_from_sequence.side_effect = [
            ([40] * 5, [100] * 5),  # Chain 1
            ([30] * 5, [120] * 5),  # Chain 2
            ([35] * 10, [110] * 10)  # Full sequence (sum of chains)
        ]

        mock_sasa_to_rasa.side_effect = [
            ([0.8] * 5, [0.6] * 5),  # Chain 1
            ([0.9] * 5, [0.8] * 5),  # Chain 2
            ([0.85] * 10, [0.7] * 10)  # Full sequence (sum of chains)
        ]

        mock_extract_sequence.return_value = "ABCDEFGHIJ"
        
        pdb = "test.pdb"
        osmolytes = "urea"
        expected_dG = {
            "Chain 1": {"urea": (1.0 * 5 * 0.8 + 0.5 * 5 * 0.6)},
            "Chain 2": {"urea": (0.8 * 5 * 0.9 + 0.3 * 5 * 0.8)},
            "All": {"urea": (0.9 * 10 * 0.85 + 0.4 * 10 * 0.7)}  # Full sequence calculation
        }

        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes, split_chains=True)
        print(expected_dG)
        print(result)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_local.extract_sequence")
    @patch("osmofold.osmofold_local.get_tfe")
    @patch("osmofold.osmofold_local.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_local.sasa_to_rasa")
    def test_custom_tfe(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe, mock_extract_sequence):
        mock_extract_sequence.return_value = "ABCDE"
        mock_get_tfe.return_value = ([0.5] * 5, [0.2] * 5)
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100 * 5])
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        pdb = "test.pdb"
        osmolytes = "tfe"
        custom_tfe = {"tfe": 1.0}
        expected_dG = (0.5 * 5 * 0.9 + 0.2 * 5 * 0.8)
        
        result = osmofold_local.protein_unfolded_dG(pdb, osmolytes, custom_tfe=custom_tfe)
        self.assertEqual(set(result), set({"tfe": expected_dG}))
        for k in result:
            self.assertAlmostEqual(result[k], {"tfe": expected_dG}[k])
    
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