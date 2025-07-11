import unittest
from unittest.mock import patch
import numpy as np
import os
from osmofold import osmofold_lite

class TestReadFasta(unittest.TestCase):
    def setUp(self):
        """Create temporary FASTA files for testing."""
        self.test_files = {
            "multiple_sequences": ">seq1\nATGC\n>seq2\nCGTA\n>seq3\nTTAA\n",
            "single_sequence": ">seq1\nATGCGTATTA\n",
            "empty_file": "",
            "irregular_formatting": ">seq1\nATGCGT\n\n\n>seq2\n\nCGTATA\n\n",
        }
        self.file_paths = {}
        for name, content in self.test_files.items():
            path = f"{name}.fasta"
            with open(path, "w") as f:
                f.write(content)
            self.file_paths[name] = path

    def tearDown(self):
        """Remove temporary FASTA files after testing."""
        for path in self.file_paths.values():
            os.remove(path)

    def test_multiple_sequences(self):
        sequences = osmofold_lite.read_fasta(self.file_paths["multiple_sequences"])
        self.assertEqual(sequences, ["ATGC", "CGTA", "TTAA"])

    def test_single_sequence(self):
        sequences = osmofold_lite.read_fasta(self.file_paths["single_sequence"])
        self.assertEqual(sequences, ["ATGCGTATTA"])

    def test_empty_file(self):
        sequences = osmofold_lite.read_fasta(self.file_paths["empty_file"])
        self.assertEqual(sequences, [])

    def test_irregular_formatting(self):
        sequences = osmofold_lite.read_fasta(self.file_paths["irregular_formatting"])
        self.assertEqual(sequences, ["ATGCGT", "CGTATA"])

    def test_nonexistent_file(self):
        with self.assertRaises(FileNotFoundError):
            osmofold_lite.read_fasta("nonexistent.fasta")

class TestProteinUnfoldedDGLite(unittest.TestCase):
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_single_osmolyte(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe):
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)  # Mocked TFE values
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100] * 5)
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        seq = "ABCDE"
        osmolytes = "urea"
        expected_dG = 1.0 * 0.9 * 5 + 0.5 * 0.8 * 5
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_multiple_osmolytes(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe):
        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5), ([0.8] * 5, [0.3] * 5)]
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100] * 5)
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        seq = "ABCDE"
        osmolytes = ["urea", "tmao"]
        expected_dG = {
            "urea": 1.0 * 0.9 * 5 + 0.5 * 0.8 * 5,
            "tmao": 0.8 * 0.9 * 5 + 0.3 * 0.8 * 5
        }
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_custom_tfe(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe):
        mock_get_tfe.return_value = ([0.5] * 5, [0.2] * 5)
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100] * 5)
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        seq = "ABCDE"
        osmolytes = "tfe"
        custom_tfe = {"tfe": 1.0}
        expected_dG = 0.5 * 0.9 * 5 + 0.2 * 0.8 * 5
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes, custom_tfe=custom_tfe)
        self.assertEqual(result, {"tfe": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_different_concentration(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe):
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100] * 5)
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        seq = "ABCDE"
        osmolytes = "urea"
        concentration = 2.0
        expected_dG = (1.0 * 0.9 * 5 + 0.5 * 0.8 * 5) * 2.0
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes, concentration=concentration)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    def test_empty_sequence(self, mock_get_tfe):
        mock_get_tfe.return_value = ([], [])
        seq = ""
        osmolytes = "urea"
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
        self.assertEqual(result, {"urea": 0})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.get_unfolded_sasa_from_sequence")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_case_insensitivity(self, mock_sasa_to_rasa, mock_get_unfolded_sasa_from_sequence, mock_get_tfe):
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        mock_get_unfolded_sasa_from_sequence.return_value = ([38] * 5, [100] * 5)
        mock_sasa_to_rasa.return_value = ([0.9] * 5, [0.8] * 5)
        seq = "ABCDE"
        osmolytes = "UREA"  # Uppercase input
        expected_dG = 1.0 * 0.9 * 5 + 0.5 * 0.8 * 5
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})  # Should be converted to lowercase
    
    @patch("osmofold.osmofold_lite.get_tfe")
    def test_invalid_osmolyte(self, mock_get_tfe):
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        seq = "ABCDE"
        osmolytes = 123  # Invalid type
        
        with self.assertRaises(TypeError):
            osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
    
    @patch("osmofold.osmofold_lite.get_tfe")
    def test_no_osmolytes(self, mock_get_tfe):
        seq = "ABCDE"
        osmolytes = []
        
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes)
        self.assertEqual(result, {})  # Expecting an empty dictionary

class TestProteinFoldedDGLite(unittest.TestCase):

    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_single_osmolyte(self, mock_sasa_to_rasa, mock_get_tfe):
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "urea"
        expected_dG = (1.0 * 0.8 + 0.5 * 0.5) * 5
        
        result = osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_multiple_osmolytes(self, mock_sasa_to_rasa, mock_get_tfe):
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.side_effect = [([1.0] * 5, [0.5] * 5), ([0.8] * 5, [0.3] * 5)]
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = ["urea", "tmao"]
        expected_dG = {
            "urea": (1.0 * 0.8 + 0.5 * 0.5) * 5,
            "tmao": (0.8 * 0.8 + 0.3 * 0.5) * 5
        }
        
        result = osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, expected_dG)
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_custom_tfe(self, mock_sasa_to_rasa, mock_get_tfe):
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.return_value = ([0.5] * 5, [0.2] * 5)
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "tfe"
        custom_tfe = {"tfe": 1.0}
        expected_dG = (0.5 * 0.8 + 0.2 * 0.5) * 5
        
        result = osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, custom_tfe=custom_tfe)
        self.assertEqual(result, {"tfe": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_different_concentration(self, mock_sasa_to_rasa, mock_get_tfe):
        mock_sasa_to_rasa.return_value = ([0.8] * 5, [0.5] * 5)
        mock_get_tfe.return_value = ([1.0] * 5, [0.5] * 5)
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "urea"
        concentration = 2.0
        expected_dG = (1.0 * 0.8 + 0.5 * 0.5) * 5 * 2.0
        
        result = osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, concentration=concentration)
        self.assertEqual(result, {"urea": expected_dG})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_empty_sequence(self, mock_sasa_to_rasa, mock_get_tfe):
        mock_sasa_to_rasa.return_value = ([], [])
        mock_get_tfe.return_value = ([], [])
        seq = ""
        backbone_sasa = []
        sidechain_sasa = []
        osmolytes = "urea"
        
        result = osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, {"urea": 0})
    
    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_no_osmolytes(self, mock_sasa_to_rasa, mock_get_tfe):
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = []
        
        with self.assertRaises(Exception) as context:
            osmofold_lite.protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(str(context.exception), "not enough values to unpack (expected 2, got 0)")

class TestProteinDdGFoldingLite(unittest.TestCase):

    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_single_osmolyte(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0}
        mock_unfolded_dG.return_value = {"urea": -5.0}
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "urea"
        expected_ddG = -10.0 - (-5.0)
        
        result = osmofold_lite.protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, {"urea": expected_ddG})
    
    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_multiple_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0, "tmao": -12.0}
        mock_unfolded_dG.return_value = {"urea": -5.0, "tmao": -6.0}
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = ["urea", "tmao"]
        expected_ddG = {"urea": -5.0, "tmao": -6.0}
        
        result = osmofold_lite.protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, expected_ddG)
    
    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_triplet_output(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"urea": -10.0}
        mock_unfolded_dG.return_value = {"urea": -5.0}
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "urea"
        expected_result = {"urea": (-10.0, -5.0, -5.0)}
        
        result = osmofold_lite.protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, triplet=True)
        self.assertEqual(result, expected_result)
    
    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_custom_tfe(self, mock_unfolded_dG, mock_folded_dG):
        mock_folded_dG.return_value = {"tfe": -8.0}
        mock_unfolded_dG.return_value = {"tfe": -4.0}
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = "tfe"
        custom_tfe = {"tfe": 1.0}
        expected_ddG = -8.0 - (-4.0)
        
        result = osmofold_lite.protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, custom_tfe=custom_tfe)
        self.assertEqual(result, {"tfe": expected_ddG})
    
    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_no_osmolytes(self, mock_unfolded_dG, mock_folded_dG):
        seq = "ABCDE"
        backbone_sasa = [1.0] * 5
        sidechain_sasa = [1.0] * 5
        osmolytes = []
        
        result = osmofold_lite.protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes)
        self.assertEqual(result, {})  # Expecting an empty dictionary

if __name__ == "__main__":
    unittest.main()