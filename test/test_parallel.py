from osmofold import parallel
import unittest
from unittest.mock import patch, MagicMock, mock_open
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import json

class TestProcessPDB(unittest.TestCase):
    
    @patch("os.path.join", return_value="mocked_path.pdb")
    @patch("osmofold.parallel.extract_sequences_by_chains", return_value=["SEQ1", "SEQ2"])
    @patch("osmofold.parallel.extract_sequence", return_value="SEQ")
    @patch("osmofold.parallel.protein_ddG_folding")
    def test_process_pdb_split_chains(
        self, mock_ddG_folding, mock_extract_seq, mock_extract_seq_chains, mock_path_join
    ):
        mock_ddG_folding.return_value = {
            "Chain 1": {"Osm1": (1.0, 2.0, 3.0)},
            "Chain 2": {"Osm1": (4.0, 5.0, 6.0)},
        }
        directory = "mock_dir"
        osmolytes = ["Osm1"]
        custom_tfe = {}
        pdb_file = "test.pdb"

        filename, result = parallel.process_pdb(directory, osmolytes, custom_tfe, pdb_file, split_chains=True)

        self.assertEqual(filename, pdb_file)
        self.assertIn("All", result)
        self.assertIn("Chain 1", result)
        self.assertIn("Chain 2", result)
        self.assertEqual(result["All"]["osmolytes"]["Osm1"], {"dG_Folded": 5.0, "dG_Unfolded": 7.0, "ddG_Folding": 9.0})
    
    @patch("os.path.join", return_value="mocked_path.pdb")
    @patch("osmofold.parallel.extract_sequence", return_value="SEQ")
    @patch("osmofold.parallel.protein_ddG_folding", return_value={"All": {"Osm1": (1.0, 2.0, 3.0)}})
    def test_process_pdb_no_split_chains(
        self, mock_ddG_folding, mock_extract_seq, mock_path_join
    ):
        directory = "mock_dir"
        osmolytes = ["Osm1"]
        custom_tfe = {}
        pdb_file = "test.pdb"

        filename, result = parallel.process_pdb(directory, osmolytes, custom_tfe, pdb_file, split_chains=False)

        self.assertEqual(filename, pdb_file)
        self.assertIn("All", result)
        self.assertEqual(result["All"]["osmolytes"]["Osm1"], {"dG_Folded": 1.0, "dG_Unfolded": 2.0, "ddG_Folding": 3.0})
    
    @patch("os.path.join", return_value="mocked_path.pdb")
    @patch("osmofold.parallel.extract_sequence", side_effect=Exception("Mocked error"))
    def test_process_pdb_error_handling(self, mock_extract_seq, mock_path_join):
        directory = "mock_dir"
        osmolytes = ["Osm1"]
        custom_tfe = {}
        pdb_file = "test.pdb"

        filename, result = parallel.process_pdb(directory, osmolytes, custom_tfe, pdb_file, split_chains=False)

        self.assertEqual(filename, pdb_file)
        self.assertIn("error", result)
        self.assertEqual(result["error"], "Mocked error")

def mock_process_pdb(directory, osmolytes, custom_tfe, pdb_file, concentration=1.0, split_chains=False):
    if pdb_file == "file1.pdb":
        return "file1.pdb", {"error": "Mocked error"}
    return "file2.pdb", {"All": {"osmolytes": {"Osm1": {"dG_Folded": 2.0}}}}

class TestBatchProcessPDBs(unittest.TestCase):

    @patch('osmofold.parallel.os.path.isdir')
    @patch('osmofold.parallel.os.listdir')
    @patch('osmofold.parallel.process_pdb')
    @patch('osmofold.parallel.save_results_to_csv')
    @patch('osmofold.parallel.ProcessPoolExecutor')
    @patch('osmofold.parallel.as_completed')
    def test_batch_process_pdbs_success(self, mock_as_completed, mock_Executor, mock_save_csv, mock_process_pdb, mock_listdir, mock_isdir):
        # Mocking the return values for each function
        mock_isdir.return_value = True  # Directory exists
        mock_listdir.return_value = ['file1.pdb', 'file2.pdb']  # Two PDB files in the directory

        # Mock the result of process_pdb
        mock_process_pdb.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Create a mock future
        mock_future = MagicMock()
        mock_future.result.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Mock the return value of as_completed to return the mock future
        mock_as_completed.return_value = [mock_future]  # Simulate a completed future

        # Mock the ProcessPoolExecutor to return a mocked executor
        mock_executor = MagicMock(spec=ProcessPoolExecutor)
        mock_Executor.return_value = mock_executor

        # Mocking the behavior of `submit`
        mock_executor.submit.return_value = mock_future
        mock_executor.__enter__.return_value = mock_executor

        # Call the function
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        results = parallel.batch_process_pdbs(directory, osmolytes, num_workers=1, save_csv=False)

        # Assert that the results are correct
        self.assertIn('file1.pdb', results)
        self.assertEqual(results['file1.pdb'], {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Check that `submit` was called for each pdb file
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file2.pdb', 1.0, False)

        # Ensure that the mock future's result method was called
        mock_future.result.assert_called_once()
    
        # Ensure `save_results_to_csv` is not called (because save_csv=False)
        mock_save_csv.assert_not_called()

    @patch('osmofold.parallel.os.path.isdir')
    @patch('osmofold.parallel.os.listdir')
    @patch('osmofold.parallel.process_pdb')
    @patch('osmofold.parallel.save_results_to_csv')
    @patch('osmofold.parallel.ProcessPoolExecutor')
    def test_batch_process_pdbs_directory_not_exist(self, mock_Executor, mock_save_csv, mock_process_pdb, mock_listdir, mock_isdir):
        # Simulate that the directory doesn't exist
        mock_isdir.return_value = False
        
        # Call the function and assert that it raises a NotADirectoryError
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        
        with self.assertRaises(NotADirectoryError):
           parallel.batch_process_pdbs(directory, osmolytes)

    @patch('osmofold.parallel.os.path.isdir')
    @patch('osmofold.parallel.os.listdir')
    @patch('osmofold.parallel.process_pdb')
    @patch('osmofold.parallel.save_results_to_csv')
    @patch('osmofold.parallel.ProcessPoolExecutor')
    def test_batch_process_pdbs_no_pdb_files(self, mock_Executor, mock_save_csv, mock_process_pdb, mock_listdir, mock_isdir):
        # Simulate that the directory contains no PDB files
        mock_isdir.return_value = True
        mock_listdir.return_value = []  # No PDB files in the directory
        
        # Call the function and assert that it raises a FileNotFoundError
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        
        with self.assertRaises(FileNotFoundError):
           parallel.batch_process_pdbs(directory, osmolytes)

    @patch('osmofold.parallel.os.path.isdir')
    @patch('osmofold.parallel.os.listdir')
    @patch('osmofold.parallel.process_pdb')
    @patch('osmofold.parallel.save_results_to_csv')
    @patch('osmofold.parallel.ProcessPoolExecutor')
    @patch('osmofold.parallel.as_completed')
    def test_batch_process_pdbs_save_csv(self, mock_as_completed, mock_Executor, mock_save_csv, mock_process_pdb, mock_listdir, mock_isdir):
        # Mocking the return values for each function
        mock_isdir.return_value = True  # Directory exists
        mock_listdir.return_value = ['file1.pdb', 'file2.pdb']  # Two PDB files in the directory

        # Mock the result of process_pdb
        mock_process_pdb.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Create a mock future
        mock_future = MagicMock()
        mock_future.result.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Mock the return value of as_completed to return the mock future
        mock_as_completed.return_value = [mock_future]  # Simulate a completed future

        # Mock the ProcessPoolExecutor to return a mocked executor
        mock_executor = MagicMock(spec=ProcessPoolExecutor)
        mock_Executor.return_value = mock_executor

        # Mocking the behavior of `submit`
        mock_executor.submit.return_value = mock_future
        mock_executor.__enter__.return_value = mock_executor

        # Call the function with save_csv=True
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        results = parallel.batch_process_pdbs(directory, osmolytes, num_workers=1, save_csv=True)

        # Assert that the results are correct
        self.assertIn('file1.pdb', results)
        self.assertEqual(results['file1.pdb'], {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Ensure that `save_results_to_csv` was called with the correct results
        mock_save_csv.assert_called_once_with(results)

        # Check that `submit` was called for each pdb file
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file2.pdb', 1.0, False)

        # Ensure that the mock future's result method was called
        mock_future.result.assert_called_once()
    
    @patch('osmofold.parallel.os.path.isdir')
    @patch('osmofold.parallel.os.listdir')
    @patch('osmofold.parallel.process_pdb')
    @patch('osmofold.parallel.save_results_to_csv')
    @patch('osmofold.parallel.ProcessPoolExecutor')
    @patch('osmofold.parallel.as_completed')
    def test_batch_process_pdbs_multiple_workers(self, mock_as_completed, mock_Executor, mock_save_csv, mock_process_pdb, mock_listdir, mock_isdir):
        # Mocking the return values for each function
        mock_isdir.return_value = True  # Directory exists
        mock_listdir.return_value = ['file1.pdb', 'file2.pdb']  # Two PDB files in the directory

        # Mock the result of process_pdb
        mock_process_pdb.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Create a mock future
        mock_future = MagicMock()
        mock_future.result.return_value = ('file1.pdb', {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Mock the return value of as_completed to return the mock future
        mock_as_completed.return_value = [mock_future]  # Simulate a completed future

        # Mock the ProcessPoolExecutor to return a mocked executor
        mock_executor = MagicMock(spec=ProcessPoolExecutor)
        mock_Executor.return_value = mock_executor

        # Mocking the behavior of `submit`
        mock_executor.submit.return_value = mock_future
        mock_executor.__enter__.return_value = mock_executor

        # Call the function with multiple workers
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        results = parallel.batch_process_pdbs(directory, osmolytes, num_workers=2, save_csv=False)

        # Assert that the results are correct
        self.assertIn('file1.pdb', results)
        self.assertEqual(results['file1.pdb'], {'All': {'osmolyte1': (1.0, 2.0, 3.0)}})

        # Check that `submit` was called for each pdb file
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, None, 'file2.pdb', 1.0, False)

        # Ensure that the mock future's result method was called
        mock_future.result.assert_called_once()

        # Ensure `save_results_to_csv` is not called (because save_csv=False)
        mock_save_csv.assert_not_called()

class TestSaveResultsToCSV(unittest.TestCase):

    @patch("builtins.open", new_callable=mock_open)
    @patch("osmofold.parallel.csv.writer")
    @patch("osmofold.parallel.datetime.datetime")
    def test_save_results_success(self, mock_datetime, mock_csv_writer, mock_file):
        # Mock datetime to ensure a predictable filename
        mock_datetime.now.return_value.strftime.return_value = "20250207_123456"
        
        results = {
            "protein1.pdb": {
                "Chain 1": {
                    "protein_length": 150,
                    "osmolytes": {
                        "Urea": {"dG_Folded": -5.0, "dG_Unfolded": -3.0, "ddG_Folding": -2.0}
                    }
                },
                "All": {
                    "protein_length": 150,
                    "osmolytes": {
                        "Urea": {"dG_Folded": -5.2, "dG_Unfolded": -3.1, "ddG_Folding": -2.1}
                    }
                }
            }
        }

        parallel.save_results_to_csv(results)

        # Check the correct filename was used
        expected_filename = "batch_results_ddg_20250207_123456.csv"
        mock_file.assert_called_once_with(expected_filename, mode="w", newline="")

        # Check if the correct rows were written
        mock_csv_writer.return_value.writerow.assert_any_call(["PDB name", "Chain", "Protein Length", "Osmolyte", "dG_Unfolded", "dG_Folded", "ddG_Folding", "Error"])
        mock_csv_writer.return_value.writerow.assert_any_call(["protein1.pdb", "Chain 1", 150, "Urea", -3.0, -5.0, -2.0, "N/A"])
        mock_csv_writer.return_value.writerow.assert_any_call(["protein1.pdb", "All", 150, "Urea", -3.1, -5.2, -2.1, "N/A"])

    @patch("builtins.open", new_callable=mock_open)
    @patch("osmofold.parallel.csv.writer")
    @patch("osmofold.parallel.datetime.datetime")
    def test_save_results_with_errors(self, mock_datetime, mock_csv_writer, mock_file):
        mock_datetime.now.return_value.strftime.return_value = "20250207_123456"
        
        results = {
            "protein2.pdb": {
                "error": "Failed to process file"
            },
            "protein3.pdb": {
                "Chain A": {
                    "error": "Chain processing error"
                }
            }
        }

        parallel.save_results_to_csv(results)

        # Verify error rows were written correctly
        mock_csv_writer.return_value.writerow.assert_any_call(["protein2.pdb", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "Failed to process file"])
        mock_csv_writer.return_value.writerow.assert_any_call(["protein3.pdb", "Chain A", "N/A", "N/A", "N/A", "N/A", "N/A", "Chain processing error"])

    @patch("builtins.open", new_callable=mock_open)
    @patch("osmofold.parallel.csv.writer")
    @patch("osmofold.parallel.datetime.datetime")
    def test_save_results_empty(self, mock_datetime, mock_csv_writer, mock_file):
        mock_datetime.now.return_value.strftime.return_value = "20250207_123456"
        
        results = {}

        parallel.save_results_to_csv(results)

        # Check that only the header is written
        mock_csv_writer.return_value.writerow.assert_called_once_with(["PDB name", "Chain", "Protein Length", "Osmolyte", "dG_Unfolded", "dG_Folded", "ddG_Folding", "Error"])

if __name__ == '__main__':
    unittest.main()
