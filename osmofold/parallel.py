import os
import json
import numpy as np
from .osmofold_local import protein_ddG_folding, extract_sequence, extract_sequences_by_chains
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime
import csv

def process_pdb(directory, osmolytes, backbone, custom_tfe, pdb_file, concentration=1.0, split_chains=False):
    pdb_path = os.path.join(directory, pdb_file)
    try:
        # Extract chain-specific or full-sequence data
        if split_chains:
            chain_sequences = extract_sequences_by_chains(pdb_path)  # List of sequences, one per chain
            chain_labels = [f"Chain {i+1}" for i in range(len(chain_sequences))]  # Assign labels like "Chain 1", "Chain 2", ...
        else:
            full_sequence = extract_sequence(pdb_path)
            chain_sequences = [full_sequence]  # Wrap in a list for uniform processing
            chain_labels = ["All"]  # Single label for the full sequence

        results = {}

        # Call protein_ddG_folding
        osmolyte_results = protein_ddG_folding(
            pdb_path,
            osmolytes=osmolytes,
            backbone=backbone,
            custom_tfe=custom_tfe,
            triplet=True,
            concentration=concentration,
            split_chains=split_chains,
        )

        # Process each chain sequence
        for label, sequence in zip(chain_labels, chain_sequences):
            chain_length = len(sequence)
            if label not in results:
                results[label] = {"protein_length": chain_length, "osmolytes": {}}

            if "All" not in results:
                results["All"] = {"protein_length": 0, "osmolytes": {}}

            if split_chains:
                results["All"]["protein_length"] += chain_length

            chain_osmolyte_results = osmolyte_results.get(label, {}) if split_chains else osmolyte_results.get('All', {})

            for osmolyte in osmolytes:
                # Ensure you're accessing osmolyte data nested under "All"
                dG_folded, dG_unfolded, ddG = chain_osmolyte_results.get(osmolyte, (None, None, None))
                results[label]["osmolytes"][osmolyte] = {
                    "dG_Folded": dG_folded,
                    "dG_Unfolded": dG_unfolded,
                    "ddG_Folding": ddG,
                }

                if osmolyte not in results["All"]["osmolytes"]:
                    results["All"]["osmolytes"][osmolyte] = {"dG_Folded": 0, "dG_Unfolded": 0, "ddG_Folding": 0}
                
                if split_chains:
                    results["All"]["osmolytes"][osmolyte]["dG_Folded"] += dG_folded
                    results["All"]["osmolytes"][osmolyte]["dG_Unfolded"] += dG_unfolded
                    results["All"]["osmolytes"][osmolyte]["ddG_Folding"] += ddG

        results = {k: results[k] for k in sorted(results.keys(), key=lambda x: (x == "All"))}
        return pdb_file, results
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return pdb_file, {"error": str(e)}

def save_results_to_csv(results):
    """
    Saves batch processing results to a CSV file.

    Parameters:
        results (dict): Batch processing results.
    """
    output_file = f"batch_results_ddg_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(["PDB name", "Chain", "protein length", "osmolyte", "dG_Unfolded", "dG_Folded", "ddG_Folding", "error"])

        for pdb_file, pdb_results in results.items():
            if "error" in pdb_results:
                writer.writerow([pdb_file, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", pdb_results["error"]])
                continue

            for chain, chain_results in pdb_results.items():
                if "error" in chain_results:
                    writer.writerow([pdb_file, chain, "N/A", "N/A", "N/A", "N/A", "N/A", chain_results["error"]])
                    continue

                protein_length = chain_results.get("protein_length", "N/A")
                for osmolyte, values in chain_results.get("osmolytes", {}).items():
                    dG_Folded = values.get("dG_Folded", "N/A")
                    dG_Unfolded = values.get("dG_Unfolded", "N/A")
                    ddG_Folding = values.get("ddG_Folding", "N/A")
                    writer.writerow([pdb_file, chain, protein_length, osmolyte, dG_Unfolded, dG_Folded, ddG_Folding, "N/A"])

    print(f"Results saved to {output_file}")

def batch_process_pdbs(directory, osmolytes, backbone=True, custom_tfe=None, concentration=1.0, save_csv=True, num_workers=1, split_chains=False):
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
        futures = {executor.submit(process_pdb, directory, osmolytes, backbone, custom_tfe, pdb_file, concentration, split_chains): pdb_file for pdb_file in pdb_files}
        for future in as_completed(futures):
            pdb_file, result = future.result()
            results[pdb_file] = result

    # Debug: Print the results for verification
    print("Final Results:\n", json.dumps(results, indent=2))

    if save_csv:
        save_results_to_csv(results)
    return results

