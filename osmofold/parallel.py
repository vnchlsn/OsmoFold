import os
import json
import numpy as np
from .osmofold_local import protein_folded_dG, protein_unfolded_dG, protein_ddG_folding
from concurrent.futures import ProcessPoolExecutor

def process_pdb(directory, analysis, osmolytes, backbone, custom_tfe, pdb_file):
    pdb_path = os.path.join(directory, pdb_file)
    try:
        if analysis == "folded":
            result = protein_folded_dG(pdb_path, osmolytes, backbone, custom_tfe)
        elif analysis == "unfolded":
            result = protein_unfolded_dG(pdb_path, osmolytes, backbone, custom_tfe)
        elif analysis == "ddg":
            result = protein_ddG_folding(pdb_path, osmolytes, backbone, custom_tfe)
        return pdb_file, result
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return pdb_file, None

def batch_process_pdbs(directory, osmolytes, analysis="ddg", backbone=True, custom_tfe=None, save_csv=False, num_workers=1):
    """
    Batch processes PDB files in a directory for specified osmolyte analyses with parallel processing.
    """
    if analysis not in ["folded", "unfolded", "ddg"]:
        raise ValueError("Invalid analysis type. Choose 'folded', 'unfolded', or 'ddg'.")

    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in directory: {directory}")

    # Pass all needed parameters to process_pdb
    results = {}
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Create futures with direct process_pdb calls
        futures = [
            executor.submit(process_pdb, directory, analysis, osmolytes, backbone, custom_tfe, pdb_file)
            for pdb_file in pdb_files
        ]
        # Gather results as they complete
        for future in futures:
            pdb_file, result = future.result()
            results[pdb_file] = result

    if save_csv:
        save_results_to_csv(results, analysis)
    return results

def save_results_to_csv(results, analysis):
    """
    Saves batch processing results to a CSV file.

    Parameters:
        results (dict): Batch processing results.
        analysis (str): Type of analysis performed (used for naming the file).
    """
    import csv
    output_file = f"batch_results_{analysis}.csv"
    
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(["PDB File", "Osmolyte", "Result"])
        
        for pdb_file, osmolyte_results in results.items():
            if osmolyte_results is None:
                writer.writerow([pdb_file, "N/A", "Error"])
                continue
            for osmolyte, value in osmolyte_results.items():
                writer.writerow([pdb_file, osmolyte, value])

    print(f"Results saved to {output_file}")

# Example usage
if __name__ == "__main__":
    directory = "./pdb_files"
    osmolytes = ["trehalose", "tmao"]
    analysis = "ddg"  # Choose "folded", "unfolded", or "ddg"

    results = batch_process_pdbs(directory, osmolytes, analysis, save_csv=True, num_workers=4)
    print(json.dumps(results, indent=2))
