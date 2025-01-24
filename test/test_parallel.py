from osmofold import parallel
import unittest
from unittest.mock import patch, MagicMock
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import json

class TestProcessPDB(unittest.TestCase):

    @patch('osmofold.parallel.extract_sequences_by_chains')
    @patch('osmofold.parallel.extract_sequence')
    @patch('osmofold.parallel.protein_ddG_folding')
    def test_process_pdb_split_chains(self, mock_ddg_folding, mock_extract_sequence, mock_extract_sequences_by_chains):
        # Mocking the return values
        mock_extract_sequences_by_chains.return_value = ['AAAA', 'BBBB']
        mock_extract_sequence.return_value = 'ABABAB'
        
        mock_ddg_folding.return_value = {
            'Chain 1': {'osmolyte1': (1.0, 2.0, 3.0), 'osmolyte2': (1.5, 2.5, 3.5)},
            'Chain 2': {'osmolyte1': (1.1, 2.1, 3.1), 'osmolyte2': (1.6, 2.6, 3.6)}
        }
        
        directory = '/mock/path'
        osmolytes = ['osmolyte1', 'osmolyte2']
        backbone = 'backbone_data'
        custom_tfe = 'custom_tfe_data'
        pdb_file = 'mock.pdb'
        concentration = 1.0
        split_chains = True
        
        # Call the function
        pdb_file, results = parallel.process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file, concentration, split_chains)
        
        # Test the results
        self.assertEqual(pdb_file, 'mock.pdb')
        self.assertIn('Chain 1', results)
        self.assertIn('Chain 2', results)
        self.assertIn('All', results)
        
        # Check chain-specific results
        self.assertEqual(results['Chain 1']['protein_length'], 4)
        self.assertEqual(results['Chain 2']['protein_length'], 4)
        
        # Check osmolyte results for each chain
        self.assertEqual(results['Chain 1']['osmolytes']['osmolyte1']['dG_Folded'], 1.0)
        self.assertEqual(results['Chain 1']['osmolytes']['osmolyte2']['dG_Folded'], 1.5)
        
        self.assertEqual(results['Chain 2']['osmolytes']['osmolyte1']['dG_Folded'], 1.1)
        self.assertEqual(results['Chain 2']['osmolytes']['osmolyte2']['dG_Folded'], 1.6)

        # Check combined "All" results
        self.assertEqual(results['All']['protein_length'], 8)  # 4 + 4 for two chains
        self.assertEqual(results['All']['osmolytes']['osmolyte1']['dG_Folded'], 2.1)
        self.assertEqual(results['All']['osmolytes']['osmolyte2']['dG_Folded'], 3.1)

    @patch('osmofold.parallel.extract_sequences_by_chains')
    @patch('osmofold.parallel.extract_sequence')
    @patch('osmofold.parallel.protein_ddG_folding')
    def test_process_pdb_no_split_chains(self, mock_ddg_folding, mock_extract_sequence, mock_extract_sequences_by_chains):
        # Mocking the return values
        mock_extract_sequences_by_chains.return_value = ['AAAA']
        mock_extract_sequence.return_value = 'ABABAB'
        
        mock_ddg_folding.return_value = {
            'All': {'osmolyte1': (1.0, 2.0, 3.0), 'osmolyte2': (1.5, 2.5, 3.5)}
        }
        
        directory = '/mock/path'
        osmolytes = ['osmolyte1', 'osmolyte2']
        backbone = 'backbone_data'
        custom_tfe = 'custom_tfe_data'
        pdb_file = 'mock.pdb'
        concentration = 1.0
        split_chains = False
        
        # Call the function
        pdb_file, results = parallel.process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file, concentration, split_chains)
        
        # Test the results
        self.assertEqual(pdb_file, 'mock.pdb')
        self.assertIn('All', results)
        
        # Check the "All" results
        self.assertEqual(results['All']['protein_length'], 6)  # Full sequence length (ABABAB)
        
        # Check osmolyte results for "All"
        self.assertEqual(results['All']['osmolytes']['osmolyte1']['dG_Folded'], 1.0)
        self.assertEqual(results['All']['osmolytes']['osmolyte2']['dG_Folded'], 1.5)

    @patch('osmofold.parallel.extract_sequences_by_chains')
    @patch('osmofold.parallel.extract_sequence')
    @patch('osmofold.parallel.protein_ddG_folding')
    def test_process_pdb_error_handling(self, mock_ddg_folding, mock_extract_sequence, mock_extract_sequences_by_chains):
        # Simulating an error in processing
        mock_extract_sequences_by_chains.side_effect = Exception("Mock error during sequence extraction")
        
        directory = '/mock/path'
        osmolytes = ['osmolyte1']
        backbone = 'backbone_data'
        custom_tfe = 'custom_tfe_data'
        pdb_file = 'mock.pdb'
        
        # Call the function
        pdb_file, results = parallel.process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file)
        
        # Test that error is captured
        self.assertEqual(pdb_file, 'mock.pdb')
        self.assertIn('error', results)
        self.assertEqual(results['error'], 'not enough values to unpack (expected 3, got 0)')

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
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file2.pdb', 1.0, False)

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
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file2.pdb', 1.0, False)

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
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file1.pdb', 1.0, False)
        mock_executor.submit.assert_any_call(mock_process_pdb, directory, osmolytes, True, None, 'file2.pdb', 1.0, False)

        # Ensure that the mock future's result method was called
        mock_future.result.assert_called_once()

        # Ensure `save_results_to_csv` is not called (because save_csv=False)
        mock_save_csv.assert_not_called()

if __name__ == '__main__':
    unittest.main()
