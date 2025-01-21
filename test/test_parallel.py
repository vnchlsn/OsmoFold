import unittest
from unittest.mock import mock_open, patch
from osmofold import parallel
import os

# Define a mock function that is picklable
def mock_process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file, concentration = 1.0):
    return pdb_file, {"protein_length": 10}

class TestProcessPDB(unittest.TestCase):
    @patch("osmofold.parallel.extract_sequence")
    @patch("osmofold.parallel.protein_ddG_folding")
    def test_process_pdb_success(self, mock_ddG_folding, mock_extract_sequence):
        directory = "/mock/directory"
        pdb_file = "test.pdb"
        osmolytes = ["TMAO", "Urea"]
        backbone = True
        custom_tfe = None
        pdb_path = os.path.join(directory, pdb_file)
        
        # Mocking extract_sequence
        mock_extract_sequence.return_value = "ACDEFGHIKL"
        
        # Mocking protein_ddG_folding
        mock_ddG_folding.return_value = {
            "TMAO": (-5.0, -3.0, 2.0),
            "Urea": (-4.5, -2.5, 2.0)
        }
        
        expected_output = (pdb_file, {
            "protein_length": 10,
            "osmolytes": {
                "TMAO": {"dG_Folded": -5.0, "dG_Unfolded": -3.0, "ddG_Folding": 2.0},
                "Urea": {"dG_Folded": -4.5, "dG_Unfolded": -2.5, "ddG_Folding": 2.0}
            }
        })
        
        result = parallel.process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file)
        self.assertEqual(result, expected_output)

class TestBatchProcessPDBs(unittest.TestCase):
    @patch("os.listdir", return_value=["file1.pdb", "file2.pdb"])  # Simulate PDB files
    @patch("os.path.isdir", return_value=True)  # Simulate valid directory
    def test_batch_process_pdbs_success(self, mock_isdir, mock_listdir):
        directory = "/mock/directory"
        osmolytes = ["TMAO"]

        expected_output = {
            "file1.pdb": {"protein_length": 10},
            "file2.pdb": {"protein_length": 10}
        }

        # Properly replace process_pdb with a *real* function (not a MagicMock)
        with patch("osmofold.parallel.process_pdb", new=mock_process_pdb):
            result = parallel.batch_process_pdbs(directory, osmolytes, save_csv=False, num_workers=2, concentration=1.0)

        self.assertEqual(result, expected_output)

        # Test save_results_to_csv separately (no multiprocessing here)
        with patch("osmofold.parallel.save_results_to_csv") as mock_save_csv:
            parallel.batch_process_pdbs(directory, osmolytes, save_csv=True, num_workers=1)
            mock_save_csv.assert_called_once()
    
    @patch("os.path.isdir", return_value=False)
    def test_batch_process_pdbs_invalid_directory(self, mock_isdir):
        with self.assertRaises(NotADirectoryError):
            parallel.batch_process_pdbs("/invalid/directory", ["TMAO"])
    
    @patch("os.listdir", return_value=[])
    @patch("os.path.isdir", return_value=True)
    def test_batch_process_pdbs_no_pdb_files(self, mock_isdir, mock_listdir):
        with self.assertRaises(FileNotFoundError):
            parallel.batch_process_pdbs("/mock/directory", ["TMAO"])

class TestSaveResultsToCSV(unittest.TestCase):
    @patch("builtins.open", new_callable=mock_open)
    @patch("datetime.datetime")
    def test_save_results_to_csv(self, mock_datetime, mock_file):
        # Mock datetime for predictable filenames
        mock_datetime.now.return_value.strftime.return_value = "20240116_123456"

        results = {
            "protein1.pdb": {
                "protein_length": 100,
                "osmolytes": {
                    "TMAO": {"dG_Unfolded": -5.2, "dG_Folded": -10.3, "ddG_Folding": 5.1}
                }
            },
            "protein2.pdb": {
                "error": "Failed to process"
            }
        }

        parallel.save_results_to_csv(results)

        # Check the filename format
        expected_filename = "batch_results_ddg_20240116_123456.csv"
        mock_file.assert_called_once_with(expected_filename, mode="w", newline="")

        # Retrieve written data
        file_handle = mock_file()
        written_content = file_handle.write.call_args_list

        # Expected CSV content
        expected_lines = [
            "PDB name,protein length,osmolyte,dG_Unfolded,dG_Folded,ddG_Folding,error\n",
            "protein1.pdb,100,TMAO,-5.2,-10.3,5.1,N/A\n",
            "protein2.pdb,Error,N/A,N/A,N/A,N/A,Failed to process\n"
        ]

        # Normalize line endings for cross-platform compatibility
        written_lines = [call[0][0].replace("\r\n", "\n") for call in written_content]
        
        self.assertEqual(written_lines, expected_lines)

if __name__ == "__main__":
    unittest.main()