.. _parallel-docs:

**parallel.py**
=====================

A summary of the functions in `parallel.py` (designed to batch TFE calculations) and how they work.

process_pdb()
-------------

Processes a PDB file to extract sequence data and compute the change in free energy (ΔΔG) upon folding in the presence of osmolytes.

**Arguments:**

- **`directory`**: A string specifying the path to the directory containing the PDB file.  

      Example: `"/path/to/pdb_files/"`

- **`osmolytes`**: A string containing a single osmolyte or a list of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`  

      Example: `["trehalose", "sucrose"]`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`  

- **`pdb_file`**: A string containing the name of the PDB file to be processed.  

      Example: `"protein.pdb"`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.  

- **`split_chains`**: OPTIONAL. A boolean indicating whether to process chains separately (`True`) or consider the full protein (`False`). Default is `False`.  

**Returns:**  
A tuple containing:  
1. A string representing the PDB filename.  
2. A dictionary where each key corresponds to either `"All"` (aggregated results) or a specific chain, and the corresponding value contains:  
   - `protein_length`: An integer representing the number of residues.  
   - `osmolytes`: A dictionary where each key is an osmolyte, and the value is another dictionary with folding free energy changes:  
   
      Example output:  
      
      .. code-block:: python
      
            {
                "All": {
                    "protein_length": 153,
                    "osmolytes": {
                        "trehalose": {"dG_Folded": -55.2, "dG_Unfolded": -32.7, "ddG_Folding": -22.5},
                        "sucrose": {"dG_Folded": -30.4, "dG_Unfolded": -16.7, "ddG_Folding": -13.7}
                    }
                },
                "Chain 1": {
                    "protein_length": 78,
                    "osmolytes": {
                        "trehalose": {"dG_Folded": -27.1, "dG_Unfolded": -16.4, "ddG_Folding": -10.7},
                        "sucrose": {"dG_Folded": -15.2, "dG_Unfolded": -8.3, "ddG_Folding": -6.9}
                    }
                },
                "Chain 2": {
                    "protein_length": 75,
                    "osmolytes": {
                        "trehalose": {"dG_Folded": -28.1, "dG_Unfolded": -16.3, "ddG_Folding": -11.8},
                        "sucrose": {"dG_Folded": -15.2, "dG_Unfolded": -8.4, "ddG_Folding": -6.8}
                    }
                }
            }

   If an error occurs, the function returns:  

      .. code-block:: python

            {"error": "Description of the error"}

---

save_results_to_csv()
----------------------

Saves batch processing results to a CSV file, including protein folding free energy changes (ΔΔG) for different osmolytes and any encountered errors.

**Arguments:**

- **`results`**: A dictionary containing processed results for multiple PDB files.  
  Expected structure:

      .. code-block:: python

            {
                "pdb_filename": {
                    "Chain 1": {
                        "protein_length": int,
                        "osmolytes": {
                            "Osmolyte1": {"dG_Folded": float, "dG_Unfolded": float, "ddG_Folding": float},
                            ...
                        }
                    },
                    "All": { ... },
                    "error": "Error message"  # (Optional, if an error occurred)
                },
                ...
            }

**Outputs:**  
- A CSV file named in the format `batch_results_ddg_YYYYMMDD_HHMMSS.csv`, stored in the current directory.  
- Prints the filename of the saved results.  

The generated CSV file contains the following columns:

      - **PDB name**: The filename of the processed PDB.  
      - **Chain**: The protein chain identifier or `"All"` for aggregated results.  
      - **Protein Length**: The number of residues in the sequence.  
      - **Osmolyte**: The name of the osmolyte analyzed.  
      - **dG_Unfolded**: The free energy of the unfolded state (in kcal/mol).  
      - **dG_Folded**: The free energy of the folded state (in kcal/mol).  
      - **ddG_Folding**: The computed free energy change upon folding (ΔΔG).  
      - **Error**: Any error encountered during processing (if applicable).  

**Example CSV Output:**  

      .. code-block:: csv

            PDB name,Chain,Protein Length,Osmolyte,dG_Unfolded,dG_Folded,ddG_Folding,Error
            protein1.pdb,Chain 1,78,trehalose,-16.4,-27.1,-10.7,N/A
            protein1.pdb,Chain 1,78,sucrose,-8.3,-15.2,-6.9,N/A
            protein1.pdb,All,153,trehalose,-32.7,-55.2,-22.5,N/A
            protein1.pdb,All,153,sucrose,-16.7,-30.4,-13.7,N/A
            protein2.pdb,N/A,N/A,N/A,N/A,N/A,N/A,"Error: File not found"

---

batch_process_pdbs()
----------------------

Batch processes PDB files in a directory to compute folding free energy changes (ΔΔG) for specified osmolytes using parallel processing.

**Arguments:**

- **`directory`**: A string specifying the path to the directory containing PDB files.  

      Example: `"/path/to/pdb_files/"`

- **`osmolytes`**: A string containing a single osmolyte or a list of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`  

      Example: `["trehalose", "sucrose"]`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`  

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.  

- **`save_csv`**: OPTIONAL. A boolean indicating whether to save results to a CSV file. Default is `True`.  

- **`num_workers`**: OPTIONAL. An integer specifying the number of parallel workers to use for processing. Default is `1` (sequential processing).  

- **`split_chains`**: OPTIONAL. A boolean indicating whether to process chains separately (`True`) or consider the full protein (`False`). Default is `False`.  

**Returns:**  
A dictionary mapping each PDB file to its processing results.  

Example output:  

      .. code-block:: python

            {
                "protein1.pdb": {
                    "All": {
                        "protein_length": 153,
                        "osmolytes": {
                            "trehalose": {"dG_Folded": -55.2, "dG_Unfolded": -32.7, "ddG_Folding": -22.5},
                            "sucrose": {"dG_Folded": -30.4, "dG_Unfolded": -16.7, "ddG_Folding": -13.7}
                        }
                    },
                    "Chain 1": {
                        "protein_length": 78,
                        "osmolytes": {
                            "trehalose": {"dG_Folded": -27.1, "dG_Unfolded": -16.4, "ddG_Folding": -10.7},
                            "sucrose": {"dG_Folded": -15.2, "dG_Unfolded": -8.3, "ddG_Folding": -6.9}
                        }
                    },
                    "Chain 2": {
                        "protein_length": 75,
                        "osmolytes": {
                            "trehalose": {"dG_Folded": -28.1, "dG_Unfolded": -16.3, "ddG_Folding": -11.8},
                            "sucrose": {"dG_Folded": -15.2, "dG_Unfolded": -8.4, "ddG_Folding": -6.8}
                        }
                    }
                },
                "protein2.pdb": {
                    "error": "File not found"
                }
            }

**Raises:**  
- **`NotADirectoryError`**: If the provided directory does not exist.  
- **`FileNotFoundError`**: If no PDB files are found in the directory.

*If any of the functions fail to work as described, please submit a GitHub issue or contact Vincent (`vnichol2@uwyo.edu`).*