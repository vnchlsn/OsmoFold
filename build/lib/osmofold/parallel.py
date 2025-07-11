import os
import json
from .osmofold_local import protein_ddG_folding, extract_sequence, extract_sequences_by_chains
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime
import csv

def process_pdb(directory, osmolytes, custom_tfe, pdb_file, concentration=1.0, split_chains=False):
    """
    Processes a PDB file to extract sequence data and compute the change in free energy (dG) upon folding
    in the presence of osmolytes.

    Parameters:
        directory (str): Path to the directory containing the PDB file.
        osmolytes (str or list): A single osmolyte or a list of osmolytes to evaluate.
        custom_tfe (dict, optional): A dictionary with custom transfer free energy (TFE) values for osmolytes.
        pdb_file (str): The name of the PDB file to process.
        concentration (float, optional): Osmolyte concentration in molar. Default = 1.0.
        split_chains (bool, optional): Whether to process chains separately or consider the full protein. Default = False.

    Returns:
        tuple: 
            - str: The PDB filename.
            - dict: A dictionary containing sequence lengths and osmolyte-induced folding free energy changes.
              If `split_chains` is True, results are stored per chain as well as an aggregate "All" entry.

              Example structure:
              {
                  "All": {
                      "protein_length": int,
                      "osmolytes": {
                          "Osmolyte1": {"dG_Folded": float, "dG_Unfolded": float, "ddG_Folding": float},
                          ...
                      }
                  },
                  "Chain 1": {
                      "protein_length": int,
                      "osmolytes": { ... }
                  },
                  ...
              }
              If an error occurs, returns: {"error": str}
    """
    pdb_path = os.path.join(directory, pdb_file)
    try:
        # Extract sequence data, either per chain or as a full sequence
        if split_chains:
            chain_sequences = extract_sequences_by_chains(pdb_path)  # List of sequences per chain
            chain_labels = [f"Chain {i+1}" for i in range(len(chain_sequences))]  # Assign labels like "Chain 1", "Chain 2", ...
        else:
            full_sequence = extract_sequence(pdb_path)
            chain_sequences = [full_sequence]  # Wrap in a list for uniform processing
            chain_labels = ["All"]  # Single label for the full sequence

        results = {}

        # Compute ddG folding values using protein_ddG_folding
        osmolyte_results = protein_ddG_folding(
            pdb_path,
            osmolytes=osmolytes,
            custom_tfe=custom_tfe,
            triplet=True,
            concentration=concentration,
            split_chains=split_chains,
        )

        # Process extracted sequence data and organize results
        for label, sequence in zip(chain_labels, chain_sequences):
            chain_length = len(sequence)  # Determine chain length
            if label not in results:
                results[label] = {"protein_length": chain_length, "osmolytes": {}}

            # Initialize an "All" entry if not already present
            if "All" not in results:
                results["All"] = {"protein_length": 0, "osmolytes": {}}

            if split_chains:
                results["All"]["protein_length"] += chain_length  # Sum chain lengths

            # Retrieve osmolyte energy values, ensuring correct label handling
            chain_osmolyte_results = osmolyte_results if not split_chains else osmolyte_results.get(label, {})

            for osmolyte in osmolytes:
                # Extract folding energy values or default to None
                dG_folded, dG_unfolded, ddG = chain_osmolyte_results.get(osmolyte, (None, None, None))
                results[label]["osmolytes"][osmolyte] = {
                    "dG_Folded": dG_folded,
                    "dG_Unfolded": dG_unfolded,
                    "ddG_Folding": ddG,
                }

                # Ensure an entry for osmolytes exists in the "All" category
                if osmolyte not in results["All"]["osmolytes"]:
                    results["All"]["osmolytes"][osmolyte] = {"dG_Folded": 0, "dG_Unfolded": 0, "ddG_Folding": 0}
                
                # If chains are split, sum values into the "All" entry
                if split_chains:
                    results["All"]["osmolytes"][osmolyte]["dG_Folded"] += dG_folded
                    results["All"]["osmolytes"][osmolyte]["dG_Unfolded"] += dG_unfolded
                    results["All"]["osmolytes"][osmolyte]["ddG_Folding"] += ddG

        # Ensure "All" appears last in the results dictionary
        results = {k: results[k] for k in sorted(results.keys(), key=lambda x: (x == "All"))}

        return pdb_file, results

    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return pdb_file, {"error": str(e)}

def save_results_to_csv(results):
    """
    Saves batch processing results to a CSV file, including protein folding free energy changes (dG) 
    for different osmolytes and any encountered errors.

    Parameters:
        results (dict): A dictionary containing processed results for multiple PDB files. 
                        Expected structure:
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

    Outputs:
        - A CSV file named in the format `batch_results_ddg_YYYYMMDD_HHMMSS.csv`, 
          stored in the current directory.
        - Prints the filename of the saved results.
    """
    # Generate a timestamped filename for the CSV output
    output_file = f"batch_results_ddg_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)

        # Write CSV header
        writer.writerow(["PDB name", "Chain", "Protein Length", "Osmolyte", "dG_Unfolded", "dG_Folded", "ddG_Folding", "Error"])

        # Iterate through the batch results dictionary
        for pdb_file, pdb_results in results.items():
            if "error" in pdb_results:
                # Handle case where an error occurred at the file level
                writer.writerow([pdb_file, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", pdb_results["error"]])
                continue

            for chain, chain_results in pdb_results.items():
                if "error" in chain_results:
                    # Handle case where an error occurred at the chain level
                    writer.writerow([pdb_file, chain, "N/A", "N/A", "N/A", "N/A", "N/A", chain_results["error"]])
                    continue

                protein_length = chain_results.get("protein_length", "N/A")  # Extract protein length

                for osmolyte, values in chain_results.get("osmolytes", {}).items():
                    # Extract ddG folding values, using "N/A" as fallback if missing
                    dG_Folded = values.get("dG_Folded", "N/A")
                    dG_Unfolded = values.get("dG_Unfolded", "N/A")
                    ddG_Folding = values.get("ddG_Folding", "N/A")

                    # Write the processed data row to the CSV
                    writer.writerow([pdb_file, chain, protein_length, osmolyte, dG_Unfolded, dG_Folded, ddG_Folding, "N/A"])

    print(f"Results saved to {output_file}")

def batch_process_pdbs(directory, osmolytes, custom_tfe=None, concentration=1.0, save_csv=True, num_workers=1, split_chains=False):
    """
    Batch processes PDB files in a directory to compute folding free energy changes (dG) 
    for specified osmolytes using parallel processing.

    Parameters:
        directory (str): Path to the directory containing PDB files.
        osmolytes (str or list): A single osmolyte or a list of osmolytes to analyze.
        custom_tfe (dict, optional): A dictionary containing custom transfer free energy (TFE) values for osmolytes.
        concentration (float, optional): Osmolyte concentration in molar. Default = 1.0.
        save_csv (bool, optional): Whether to save results to a CSV file. Default = True.
        num_workers (int, optional): Number of parallel workers for processing. Default = 1 (sequential processing).
        split_chains (bool, optional): Whether to process chains separately. Default = False.

    Returns:
        dict: A dictionary mapping each PDB file to its processing results.
              Example structure:
              {
                  "protein1.pdb": {
                      "All": {
                          "protein_length": int,
                          "osmolytes": {
                              "Osmolyte1": {"dG_Folded": float, "dG_Unfolded": float, "ddG_Folding": float},
                              ...
                          }
                      },
                      "Chain 1": { ... },
                      "error": "Error message"  # (Optional, if an error occurred)
                  },
                  ...
              }

    Raises:
        NotADirectoryError: If the provided directory does not exist.
        FileNotFoundError: If no PDB files are found in the directory.
    """
    # Verify that the directory exists
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"The provided directory {directory} does not exist.")

    # Collect all PDB files in the directory
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in directory: {directory}")

    results = {}

    # Process PDB files in parallel using a process pool
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit processing tasks for each PDB file
        futures = {executor.submit(process_pdb, directory, osmolytes, custom_tfe, pdb_file, concentration, split_chains): pdb_file for pdb_file in pdb_files}

        # Collect results as they complete
        for future in as_completed(futures):
            pdb_file, result = future.result()
            results[pdb_file] = result

    # Debugging: Print the results for verification
    print("Final Results:\n", json.dumps(results, indent=2))

    # Save results to CSV if requested
    if save_csv:
        save_results_to_csv(results)

    return results

