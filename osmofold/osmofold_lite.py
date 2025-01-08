import numpy as np
from .osmofold_local import sasa_to_rasa, get_tfe
    
def read_fasta(fasta_path):
    """
    Reads a FASTA file and converts it into a list of sequences.

    Parameters:
        fasta_path (str): Path to the FASTA file.

    Returns:
        list: A list of sequences (as strings) from the FASTA file.
    """
    sequences = []
    with open(fasta_path, 'r') as fasta_file:
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # Save the previous sequence if any
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                # Append to the current sequence
                sequence += line
        # Add the last sequence to the list
        if sequence:
            sequences.append(sequence)
    return sequences
    
def protein_unfolded_dG_lite(seq, osmolytes, backbone=True, custom_tfe=None, concentration=1.0):
    """
    Computes the total free energy (dG) for the unfolded protein for one or multiple osmolytes.

    Parameters:
        seq (str): The sequence of interest, as one-letter code.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        backbone (bool): Whether to account for the protein backbone. Default = True.
        custom_tfe (dict, optional): A dictionary with custom TFE values for osmolytes.
        concentration (float): The concentration in molar to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is the total free energy.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    results = {}

    for osmo in osmolytes:
        osmo = osmo.lower()
        if backbone:
            osmoBack = osmo + "Back"
        else:
            osmoBack = osmo
        results[osmo] = concentration * np.sum(get_tfe(seq, osmoBack, custom_tfe))

    return results

def protein_folded_dG_lite(seq, sasa_list, osmolytes, backbone=True, custom_tfe=None, concentration=1.0):
    """
    Computes the total free energy (dG) for the unfolded protein for one or multiple osmolytes.

    Parameters:
        seq (str): The sequence of interest, as one-letter code.
        sasa_list (list): A list of the sasa values for each residue in the protein of interest.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        backbone (bool): Whether to account for the protein backbone. Default = True.
        custom_tfe (dict, optional): A dictionary with custom TFE values for osmolytes.
        concentration (float): The concentration in molar to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is the total free energy.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    results = {}
    rasa = sasa_to_rasa(seq, sasa_list)

    for osmo in osmolytes:
        osmo = osmo.lower()
        if backbone:
            osmoBack = osmo + "Back"
        else:
            osmoBack = osmo
        results[osmo] = concentration * np.sum(np.fromiter((tfe * rsa for tfe, rsa in zip(get_tfe(seq, osmoBack, custom_tfe), rasa)), dtype=float))

    return results

def protein_ddG_folding_lite(seq, sasa_list, osmolytes, backbone=True, triplet=False, custom_tfe=None, concentration=1.0):
    """
    Computes the change in free energy (dG) upon protein folding for one or multiple osmolytes.

    Parameters:
        seq (str): The sequence of interest, as one-letter code.
        sasa_list (list): A list of the sasa values for each residue in the protein of interest.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        backbone (bool): Whether to account for the protein backbone. Default = True.
        triplet (bool): Whether to return the triplet (folded, unfolded, and their difference). Default = False.
        custom_tfe (dict, optional): A dictionary with custom TFE values for osmolytes.
        concentration (float): The concentration in molar to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is either a tuple (folded_dG, unfolded_dG, dG_change) or the free energy difference.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    folded_dG = protein_folded_dG_lite(seq, sasa_list, osmolytes, backbone, custom_tfe, concentration)
    unfolded_dG = protein_unfolded_dG_lite(seq, osmolytes, backbone, custom_tfe, concentration)

    results = {}
    for osmo in osmolytes:
        if triplet:
            results[osmo] = (folded_dG[osmo], unfolded_dG[osmo], folded_dG[osmo] - unfolded_dG[osmo])
        else:
            results[osmo] = folded_dG[osmo] - unfolded_dG[osmo]

    return results