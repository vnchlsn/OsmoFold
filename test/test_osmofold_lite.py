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
    @patch('osmofold.osmofold_lite.get_tfe')
    def test_custom_tfe(self, mock_get_tfe):
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]  # Custom TFE should not matter here
        
        custom_tfe = {"tmaoBack": 1.0}
        result = osmofold_lite.protein_unfolded_dG_lite("ACD", "tmao", custom_tfe=custom_tfe)
        expected = {"tmao": 6.0}
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.get_tfe")
    def test_multiple_osmolytes_without_backbone(self, mock_get_tfe):
        seq = "ACD"
        osmolytes = ["osmolyte1", "osmolyte2"]

        # Mock `get_tfe` to return valid numeric values
        def mock_get_tfe_function(seq, osmolyte, custom_tfe=None):
            if osmolyte == "osmolyte1":
                return [0.5, 1.0, 1.5]  # Mocked TFE values for osmolyte1
            elif osmolyte == "osmolyte2":
                return [1.0, 1.0, 1.0]  # Mocked TFE values for osmolyte2
            return [0.0]  # Default return value for other cases

        mock_get_tfe.side_effect = mock_get_tfe_function

        # Call the function to test
        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolytes, backbone=False)

        # Define the expected result
        expected = {
            "osmolyte1": 3.0,  # sum([0.5, 1.0, 1.5])
            "osmolyte2": 3.0,  # sum([1.0, 1.0, 1.0])
        }

        # Assert that the result matches the expectation
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.get_tfe")
    def test_different_concentrations(self, mock_get_tfe):
        seq = "ACDEFGHIK"
        osmolyte = "osmolyte1"
        concentration = 2.0
        mock_get_tfe.return_value = [2.0, 3.0, 1.0]

        result =  osmofold_lite.protein_unfolded_dG_lite(seq, osmolyte, concentration=concentration)
        expected = {"osmolyte1": 12.0}  # 2.0 * sum([1.0, 2.0, 3.0])
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.get_tfe")
    def test_empty_sequence(self, mock_get_tfe):
        seq = ""
        osmolyte = "osmolyte1"
        mock_get_tfe.return_value = []

        result = osmofold_lite.protein_unfolded_dG_lite(seq, osmolyte)
        expected = {"osmolyte1": 0.0}
        self.assertEqual(result, expected)

    def test_unknown_osmolyte(self):
        seq = "ACDEFGHIK"
        osmolyte = "unknown"

        with self.assertRaises(TypeError):
            osmofold_lite.protein_unfolded_dG_lite(seq, osmolyte)

class TestProteinFoldedDGLite(unittest.TestCase):

    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_single_osmolyte_with_backbone(self, mock_sasa_to_rasa, mock_get_tfe):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1"]
        concentration = 1.0
        
        # Setting up mock returns
        mock_get_tfe.return_value = [1.0, 1.2, 1.3, 0.9, 1.1, 1.0, 1.3, 1.2, 1.5]  # mock TFE values
        mock_sasa_to_rasa.return_value = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2]  # mock RSA values

        # Calling the function under test
        result = osmofold_lite.protein_folded_dG_lite(seq, sasa_list, osmolytes, backbone=True)

        # Expected sum of (TFE * RSA) for osmolyte1 with the given sequence and SASA values
        expected = {
            "osmolyte1": concentration * np.sum([
                tfe * rsa
                for tfe, rsa in zip([1.0, 1.2, 1.3, 0.9, 1.1, 1.0, 1.3, 1.2, 1.5], 
                                     [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2])
            ])
        }

        # Validate the result
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_multiple_osmolytes_without_backbone(self, mock_sasa_to_rasa, mock_get_tfe):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1", "osmolyte2"]

        # Setting up mock returns for multiple osmolytes
        mock_get_tfe.side_effect = [
            [1.0, 1.1, 1.2, 1.0, 1.1, 1.0, 1.2, 1.3, 1.4],  # osmolyte1 TFE values
            [0.9, 1.0, 1.1, 1.0, 0.9, 1.1, 1.0, 0.9, 1.0]   # osmolyte2 TFE values
        ]
        mock_sasa_to_rasa.return_value = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2]  # mock RSA values

        # Calling the function under test
        result = osmofold_lite.protein_folded_dG_lite(seq, sasa_list, osmolytes, backbone=False)

        # Expected sums for osmolyte1 and osmolyte2
        expected = {
            "osmolyte1": np.sum([
                tfe * rsa
                for tfe, rsa in zip([1.0, 1.1, 1.2, 1.0, 1.1, 1.0, 1.2, 1.3, 1.4],
                                     [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2])
            ]),
            "osmolyte2": np.sum([
                tfe * rsa
                for tfe, rsa in zip([0.9, 1.0, 1.1, 1.0, 0.9, 1.1, 1.0, 0.9, 1.0],
                                     [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2])
            ])
        }

        # Validate the result
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.get_tfe")
    @patch("osmofold.osmofold_lite.sasa_to_rasa")
    def test_custom_tfe(self, mock_sasa_to_rasa, mock_get_tfe):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1"]

        # Custom TFE values for testing
        custom_tfe = {"osmolyte1Back": [2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0]}

        # Setting up mock returns
        mock_get_tfe.return_value = custom_tfe["osmolyte1Back"]
        mock_sasa_to_rasa.return_value = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2]  # mock RSA values

        # Calling the function under test
        result = osmofold_lite.protein_folded_dG_lite(seq, sasa_list, osmolytes, custom_tfe=custom_tfe, backbone=True)

        # Expected sum of (custom TFE * RSA)
        expected = {
            "osmolyte1": np.sum([
                tfe * rsa
                for tfe, rsa in zip([2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0],
                                     [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.4, 0.3, 0.2])
            ])
        }

        # Validate the result
        self.assertEqual(result, expected)

class TestProteinDdGFoldingLite(unittest.TestCase):

    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_single_osmolyte_with_backbone(self, mock_protein_unfolded_dG_lite, mock_protein_folded_dG_lite):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1"]

        # Setting up mock return values for folded and unfolded free energies
        mock_protein_folded_dG_lite.return_value = {"osmolyte1": 5.0}
        mock_protein_unfolded_dG_lite.return_value = {"osmolyte1": 3.0}

        # Calling the function under test
        result = osmofold_lite.protein_ddG_folding_lite(seq, sasa_list, osmolytes, backbone=True)

        # Expected result: folded_dG - unfolded_dG
        expected = {
            "osmolyte1": 5.0 - 3.0
        }

        # Validate the result
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_multiple_osmolytes_without_backbone(self, mock_protein_unfolded_dG_lite, mock_protein_folded_dG_lite):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1", "osmolyte2"]

        # Setting up mock return values for folded and unfolded free energies
        mock_protein_folded_dG_lite.return_value = {"osmolyte1": 5.0, "osmolyte2": 6.0}
        mock_protein_unfolded_dG_lite.return_value = {"osmolyte1": 3.0, "osmolyte2": 4.0}

        # Calling the function under test
        result = osmofold_lite.protein_ddG_folding_lite(seq, sasa_list, osmolytes, backbone=False)

        # Expected results: folded_dG - unfolded_dG for each osmolyte
        expected = {
            "osmolyte1": 5.0 - 3.0,
            "osmolyte2": 6.0 - 4.0
        }

        # Validate the result
        self.assertEqual(result, expected)

    @patch("osmofold.osmofold_lite.protein_folded_dG_lite")
    @patch("osmofold.osmofold_lite.protein_unfolded_dG_lite")
    def test_triplet_option(self, mock_protein_unfolded_dG_lite, mock_protein_folded_dG_lite):
        seq = "ACDEFGHIK"
        sasa_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
        osmolytes = ["osmolyte1"]

        # Setting up mock return values for folded and unfolded free energies
        mock_protein_folded_dG_lite.return_value = {"osmolyte1": 5.0}
        mock_protein_unfolded_dG_lite.return_value = {"osmolyte1": 3.0}

        # Calling the function under test with triplet=True
        result = osmofold_lite.protein_ddG_folding_lite(seq, sasa_list, osmolytes, triplet=True)

        # Expected result: (folded_dG, unfolded_dG, dG_change) for osmolyte1
        expected = {
            "osmolyte1": (5.0, 3.0, 5.0 - 3.0)
        }

        # Validate the result
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()