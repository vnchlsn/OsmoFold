import numpy as np
from .osmofold_local import get_unfolded_sasa_from_sequence, sasa_to_rasa, get_tfe, clean_dict
    
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
    
def protein_unfolded_dG_lite(seq, osmolytes, custom_tfe=None, concentration=1.0):
    """
    Computes the total free energy (dG) for the unfolded state of a protein in the presence of one or multiple osmolytes.

    Parameters:
        seq (str): The amino acid sequence of the protein, using one-letter codes.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        custom_tfe (dict, optional): A dictionary with custom transfer free energy (TFE) values for osmolytes.
        concentration (float): The osmolyte concentration in molar units, used to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is the total free energy of the unfolded protein.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    results = {}

    for osmo in osmolytes:
        osmo = osmo.lower()

        backbone_sasa, sidechain_sasa = get_unfolded_sasa_from_sequence(seq)
        backbone_rasa, sidechain_rasa = sasa_to_rasa(seq, backbone_sasa, sidechain_sasa)

        backbone_tfe, sidechain_tfe = get_tfe(seq, osmo, custom_tfe)

        backbone_result = concentration * np.sum(np.array([
            backbone_rasa[i] * backbone_tfe[i] for i in range(len(seq))]))
        sidechain_result = concentration * np.sum(np.array([
            sidechain_rasa[i] * sidechain_tfe[i] for i in range(len(seq))]))

        results[osmo] = backbone_result + sidechain_result

    return clean_dict(results)

def protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, custom_tfe=None, concentration=1.0):
    """
    Computes the total free energy (dG) for the folded state of a protein in the presence of one or multiple osmolytes.

    Parameters:
        seq (str): The amino acid sequence of the protein, using one-letter codes.
        backbone_sasa (list): A list of solvent-accessible surface area (SASA) values for the backbone of each residue.
        sidechain_sasa (list): A list of SASA values for the sidechain of each residue.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        custom_tfe (dict, optional): A dictionary with custom transfer free energy (TFE) values for osmolytes.
        concentration (float): The osmolyte concentration in molar units, used to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is the total free energy of the folded protein.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    results = {}
    backbone_rasa, sidechain_rasa = sasa_to_rasa(seq, backbone_sasa, sidechain_sasa)

    for osmo in osmolytes:
        osmo = osmo.lower()

        backbone_tfe, sidechain_tfe = get_tfe(seq, osmo, custom_tfe)

        backbone_result = concentration * np.sum(np.array([
            backbone_rasa[i] * backbone_tfe[i] for i in range(len(seq))]))
        sidechain_result = concentration * np.sum(np.array([
            sidechain_rasa[i] * sidechain_tfe[i] for i in range(len(seq))]))

        results[osmo] = backbone_result + sidechain_result

    return clean_dict(results)

def protein_ddG_folding_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, triplet=False, custom_tfe=None, concentration=1.0):
    """
    Computes the change in free energy (ΔΔG) of a protein conformational change in the presence of one or multiple osmolytes.

    Parameters:
        seq (str): The amino acid sequence of the protein, using one-letter codes.
        backbone_sasa (list): A list of solvent-accessible surface area (SASA) values for the backbone of each residue.
        sidechain_sasa (list): A list of SASA values for the sidechain of each residue.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        triplet (bool): Whether to return a tuple with the folded and unfolded free energies along with their difference. Default = False.
        custom_tfe (dict, optional): A dictionary with custom transfer free energy (TFE) values for osmolytes.
        concentration (float): The osmolyte concentration in molar units, used to scale the result. Default = 1.0.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is either a tuple (folded_dG, unfolded_dG, ΔG_change) if triplet=True,
              or the free energy difference (ΔG_change) if triplet=False.
    """
    if isinstance(osmolytes, str):
        osmolytes = [osmolytes]

    folded_dG = protein_folded_dG_lite(seq, backbone_sasa, sidechain_sasa, osmolytes, custom_tfe, concentration)
    unfolded_dG = protein_unfolded_dG_lite(seq, osmolytes, custom_tfe, concentration)

    results = {}
    for osmo in osmolytes:
        if triplet:
            results[osmo] = (folded_dG[osmo], unfolded_dG[osmo], folded_dG[osmo] - unfolded_dG[osmo])
        else:
            results[osmo] = folded_dG[osmo] - unfolded_dG[osmo]

    return clean_dict(results)