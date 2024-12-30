import os
import json
import numpy as np
from .osmofold_local import protein_ddG_folding, extract_sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime
import csv

def process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file):
    pdb_path = os.path.join(directory, pdb_file)
    try:
        sequence = extract_sequence(pdb_path)
        protein_length = len(sequence)
        results = {"protein_length": protein_length, "osmolytes": {}}

        osmolyte_results = protein_ddG_folding(pdb_path, osmolytes=osmolytes, backbone=backbone, custom_tfe=custom_tfe, triplet=True)
        for osmolyte in osmolytes:
            dG_unfolded, dG_folded, ddG = osmolyte_results[osmolyte]
            results["osmolytes"][osmolyte] = {
                "dG_Folded": dG_folded,
                "dG_Unfolded": dG_unfolded,
                "ddG_Folding": ddG
            }
        return pdb_file, results
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return pdb_file, {"error": str(e)}

def batch_process_pdbs(directory, osmolytes, backbone=True, custom_tfe=None, save_csv=True, num_workers=1):
    """
    Batch processes PDB files in a directory for specified osmolyte analyses with parallel processing.
    """
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"The provided directory {directory} does not exist.")

    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in directory: {directory}")

    results = {}
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_pdb, directory, osmolytes, backbone, custom_tfe, pdb_file): pdb_file for pdb_file in pdb_files}
        for future in as_completed(futures):
            pdb_file, result = future.result()
            results[pdb_file] = result

    if save_csv:
        save_results_to_csv(results)
    return results

def save_results_to_csv(results):
    """
    Saves batch processing results to a CSV file with extended details.

    Parameters:
        results (dict): Batch processing results.
    """
    output_file = f"batch_results_ddg_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(["PDB name", "protein length", "osmolyte", "dG_Unfolded", "dG_Folded", "ddG_Folding", "error"])

        for pdb_file, osmolyte_results in results.items():
            if "error" in osmolyte_results:
                # Write an error row if processing failed
                writer.writerow([pdb_file, "Error", "N/A", "N/A", "N/A", "N/A", osmolyte_results["error"]])
                continue

            protein_length = osmolyte_results.get("protein_length", "N/A")

            for osmolyte, values in osmolyte_results.get("osmolytes", {}).items():
                # Extract dG_Unfolded, dG_Folded, and ddG_Folding from the results
                dG_Unfolded = values.get("dG_Unfolded", "N/A")
                dG_Folded = values.get("dG_Folded", "N/A")
                ddG_Folding = values.get("ddG_Folding", "N/A")

                writer.writerow([pdb_file, protein_length, osmolyte, dG_Unfolded, dG_Folded, ddG_Folding, "N/A"])

    print(f"Results saved to {output_file}")

# Example usage
if __name__ == "__main__":
    directory = "./pdb_files"
    osmolytes = ["trehalose", "tmao"]

    results = batch_process_pdbs(directory, osmolytes, save_csv=True, num_workers=16)
    print(json.dumps(results, indent=2))

