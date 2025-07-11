import numpy as np
import mdtraj as md

def clean_dict(d):
    """
    Recursively converts all numpy float64 objects in a dictionary to Python floats.

    Args:
        d (dict, list, or tuple): A nested dictionary, list, or tuple that may contain numpy float64 objects.

    Returns:
        dict, list, or tuple: The same structure with numpy float64 objects converted to Python floats.
    """
    if isinstance(d, dict):
        return {k: clean_dict(v) for k, v in d.items()}
    elif isinstance(d, list):
        return [clean_dict(v) for v in d]
    elif isinstance(d, tuple):
        return tuple(clean_dict(v) for v in d)
    elif isinstance(d, np.float64):
        return float(d)
    else:
        return d

def get_max_sasa_list():
    """
    Returns a list of maximum SASA (Solvent Accessible Surface Area) values for each amino acid in a Gly-X-Gly tripeptide. Values from Miller et al. (1987).

    Returns:
        dictionary: A dictionary of values for the maximum solvent exposure of the protein backbone and sidechains in Gly-X-Gly tripeptide.
    """

    asa_max = {
        "backbone": {
            "A": 46, "R": 45, "N": 45, "D": 45, "C": 36,
            "Q": 45, "E": 45, "G": 85, "H": 43, "I": 42,
            "L": 43, "K": 44, "M": 44, "F": 43, "P": 38,
            "S": 42, "T": 44, "W": 42, "Y": 42, "V": 43
        },
        "sidechain": {
            "A": 67, "R": 196, "N": 113, "D": 106, "C": 104,
            "Q": 144, "E": 138, "G": 1.00, "H": 151, "I": 140,
            "L": 137, "K": 167, "M": 160, "F": 175, "P": 105,
            "S": 80, "T": 102, "W": 217, "Y": 187, "V": 117
        }}

    return asa_max

def get_unfolded_sasa_list():
    """
    Returns a list of approximate SASA (Solvent Accessible Surface Area) values for each amino acid in a disordered chain.

    Returns:
        dictionary: A dictionary of values for the solvent exposure of the protein backbone and sidechains in a random disordered protein.
    """

    asa_means = {
        "backbone": {
            "A": 27.85, "R": 25.05, "N": 25.15, "D": 26.00, "C": 26.35,
            "Q": 25.30, "E": 25.70, "G": 65.15, "H": 24.15, "I": 19.95,
            "L": 22.70, "K": 26.05, "M": 25.25, "F": 24.30, "P": 22.50,
            "S": 29.40, "T": 24.05, "W": 23.55, "Y": 25.60, "V": 20.40
        },
        "sidechain": {
            "A": 55.10, "R": 171.10, "N": 90.05, "D": 87.00, "C": 72.95,
            "Q": 116.85, "E": 113.35, "G": 1.00, "H": 111.50, "I": 117.10,
            "L": 109.55, "K": 150.65, "M": 122.40, "F": 129.25, "P": 87.00,
            "S": 66.50, "T": 84.25, "W": 156.55, "Y": 141.65, "V": 96.35
        }}

    return asa_means

def get_unfolded_sasa_from_sequence(sequence):
    """
    Given a single-letter amino acid sequence, returns two lists:
    - Backbone SASA values
    - Sidechain SASA values
    based on the unfolded/disordered reference state.

    Parameters:
        sequence (str): A string of single-letter amino acid codes.

    Returns:
        tuple: (backbone_sasa_list, sidechain_sasa_list)
    """
    sequence = sequence.upper()
    asa_means = get_unfolded_sasa_list()

    backbone_sasa = []
    sidechain_sasa = []

    for aa in sequence:
        if aa not in asa_means['backbone'] or aa not in asa_means['sidechain']:
            raise ValueError(f"Invalid or unsupported amino acid: '{aa}'")

        backbone_sasa.append(asa_means['backbone'][aa])
        sidechain_sasa.append(asa_means['sidechain'][aa])

    return backbone_sasa, sidechain_sasa

def amino_to_energy(amino, cosolute):
    """
    Determines the energy contribution of an amino acid based on the cosolute environment.
    
    If the cosolute has "Back" in its name, the function returns a dictionary where all amino acids map to the constant value.
    
    Parameters:
        amino (str): The amino acid one-letter code.
        cosolute (str): The cosolute environment.

    Returns:
        float or dict: The energy value for the given amino acid and cosolute, or a dictionary if "Back" is present.
    """
    cosolute_dicts = {
        "trehalose_hong": (
            {
                "A": 59.3, "F": 57.1, "L": 190.8, "I": 221.3, "V": 153.0, "P": 135.7, "M": 246.56,
                "G": 0, "S": -1.58, "T": 13.20, "Q": -23.98, "N": -123.1, "D": -67.6, "E": -86.7,
                "H": -68.0, "K": 204.3, "R": 193.9, "C": 0, "W": 35.9, "Y": -96.1
            },
            35
        ),
        "tmao": (
            {
                "A": -14.64, "F": -9.32, "L": 11.62, "I": -25.43, "V": -1.02, "P": -137.73, "M": -7.65,
                "W": -152.87, "G": 0, "S": -39.05, "T": 3.57, "Y": -114.32, "Q": 41.41, "N": 55.69,
                "D": -66.67, "E": -83.25, "H": 42.07, "K": -110.23, "R": -109.27, "C": 0
            },
            90
        ),
        "sarcosine": (
            {
                "A": 10.91, "F": -12.64, "L": 38.33, "I": 39.98, "V": 29.32, "P": -32.23, "M": 8.18,
                "W": -113.03, "G": 0, "S": -27.98, "T": -7.54, "Y": -26.37, "Q": -10.19, "N": -40.93,
                "D": -14.20, "E": -12.61, "H": -20.80, "K": -27.42, "R": -32.24, "C": 0
            },
            52
        ),
        "betaine": (
            {
                "A": 4.77, "F": -112.93, "L": -17.73, "I": -1.27, "V": -19.63, "P": -125.16, "M": -14.16, 
                "W": -369.93, "G": 0, "S": -241.85, "T": 0.33, "Y": -213.09, "Q": 7.57, "N": 33.17, 
                "D": -116.56, "E": -112.08, "H": -35.97, "K": -171.99, "R": -109.45, "C": 0
            },
            67
        ),
        "sorbitol": (
            {
                 "A": 16.57, "F": 26.38, "L": 39.07, "I": 36.90, "V": 24.65, "P": -4.48, "M": 20.97, 
                 "W": -67.23, "G": 0, "S": -1.58, "T": 13.20, "Y": -53.50, "Q": -23.98, "N": -21.21, 
                 "D": -83.88, "E": -70.05, "H": -42.45, "K": -32.47, "R": -24.65, "C": 0
            },
            35
        ),
        "sucrose": (
            {
                "A": 22.05, "F": -96.35, "L": 31.11, "I": 28.12, "V": 33.92, "P": -73.02, "M": -6.66, 
                "G": 0, "S": -2.79, "T": 20.82, "Q": -40.87, "N": -28.28, "D": -37.17, "E": -41.65, 
                "H": -118.66, "K": -39.60, "R": -79.32, "C": 0, "W": -215.27, "Y": -78.41
            },
            62
        ),
        "urea": (
            {
                "A": -4.69, "F": -83.11, "L": -54.57, "I": -38.43, "V": -21.65, "P": -17.65, "M": -48.34, 
                "G": 0, "S": -20.56, "T": 13.20, "Q": -54.81, "N": -38.79, "D": 3.55, "E": 0.62, 
                "H": -50.51, "K": -22.79, "R": -21.17, "C": 0, "W": -141.46, "Y": -45.08
            },
            -39
        ),
        "proline": (
            {
                "A": -0.07, "F": -71.26, "L": 4.77, "I": -2.72, "V": 7.86, "P": -63.96, "M": -35.12, 
                "G": 0, "S": -33.49, "T": -18.33, "Q": -32.36, "N": -17.71, "D": -90.51, "E": -89.17, 
                "H": -45.10, "K": -59.87, "R": -60.18, "C": 0, "W": -198.37, "Y": -138.41
            },
            48
        ),
        "glycerol": (
            {
                "A": 7.77, "F": 59.78, "L": -34.41, "I": 36.23, "V": -1.37, "P": -60.54, "M": 13.88, 
                "W": -122.65, "G": 0, "S": 6.31, "T": 17.54, "Y": -149.50, "Q": -2.75, "N": 51.57, 
                "D": -85.45, "E": -74.20, "H": -17.16, "K": -34, "R": -30.74, "C": 0
            },
            14
        ),
        "trehalose": (
            {
                "A": 33.25, "F": -17.88, "L": 96.17, "I": 79.66, "V": 96.79, "P": -94.67, "M": 29.19, 
                "W": -206.30, "G": 0, "S": -0.98, "T": 26.32, "Y": -80.32, "Q": -36.34, "N": 48.67, 
                "D": -96.54, "E": -85.92, "H": -98.75, "K": -50.08, "R": -50.33, "C": 0
            },
            62
        )
    }
    
    base_cosolute = cosolute.replace("Back", "").strip()
    if base_cosolute in cosolute_dicts:
        energy_dict, back_value = cosolute_dicts[base_cosolute]
        
        if "Back" in cosolute:
            return back_value  # Return a single value, not a dictionary
        
        return energy_dict.get(amino, None)
    
    return None

def extract_sequence(pdb_file):
    """
    Extracts the amino acid sequence from a PDB file.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        str: The amino acid sequence as a string.
    """
    sequence = []
    seen_residues = set()  # To track and avoid duplicates

    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                res_id = (line[21], line[22:26].strip())  # Chain ID and residue number
                if res_id not in seen_residues:
                    seen_residues.add(res_id)
                    residue_name = line[17:20].strip()
                    aa = three_to_one(residue_name)
                    sequence.append(aa)
                    
    return ''.join(sequence)

def extract_sequences_by_chains(pdb_file):
    """
    Extracts the amino acid sequences by chains from a PDB file.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        list: A list of amino acid sequences, one for each chain.
    """
    sequences = {}
    seen_residues = set()  # To track and avoid duplicates within each chain

    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]  # Chain ID
                res_id = (chain_id, line[22:26].strip())  # Chain ID and residue number
                if res_id not in seen_residues:
                    seen_residues.add(res_id)
                    residue_name = line[17:20].strip()
                    aa = three_to_one(residue_name)

                    if chain_id not in sequences:
                        sequences[chain_id] = []  # Initialize the sequence list for this chain
                    sequences[chain_id].append(aa)

    # Combine the sequences for each chain into strings
    return [''.join(seq) for seq in sequences.values()]

def three_to_one(residue):
    conversion_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 
        'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 
        'TYR': 'Y', 'VAL': 'V'
    }
    return conversion_dict.get(residue, 'X')  # 'X' for unknown residues

def get_tfe(seq, osmo, custom_tfe=None):
    """
    Computes the Total Free Energy (TFE) for a given sequence based on the osmolyte environment,
    returning both sidechain and backbone TFEs.

    Parameters:
        seq (str): The amino acid sequence.
        osmo (str): The osmolyte of interest.
        custom_tfe (dict, optional): A dictionary with custom TFE values structured as:
            {
                "backbone": { 'A': X, 'F': Y, 'L': Z, ... },
                "sidechain": { 'A': A, 'F': B, 'L': C, ... }
            }
            WILL BE USED IN PLACE OF THE OSMOLYTE LIST

    Returns:
        tuple: Two lists of TFE values for each amino acid in the sequence (sidechain, backbone).
    """
    amino_acids = 'AFLIVPMGSTQNDEHKRCWY'
    
    # Handle custom TFE if provided
    if custom_tfe:
        sidechain_tfe_dict = custom_tfe["sidechain"]
        backbone_tfe_dict = custom_tfe["backbone"]
    else:
        sidechain_tfe_dict = {aa: amino_to_energy(aa, osmo) for aa in amino_acids}
        backbone_tfe_dict = {aa: amino_to_energy(aa, osmo + "Back") for aa in amino_acids}
    
    # Compute TFE lists
    backbone_tfe = [backbone_tfe_dict.get(aa, 0) for aa in seq]
    sidechain_tfe = [sidechain_tfe_dict.get(aa, 0) for aa in seq]
    
    return backbone_tfe, sidechain_tfe

def get_pdb_info(pdb):
    """
    Loads a PDB file and calculates per-residue SASA values for backbone and sidechain atoms,
    along with the amino acid sequence.

    Parameters:
        pdb (str): Path to the PDB file.

    Returns:
        tuple of lists:
            - seq: List of single-letter amino acid codes.
            - backbone_sasa: List of backbone SASA values for each residue.
            - sidechain_sasa: List of sidechain SASA values for each residue.
    """
    traj = md.load_pdb(pdb)

    sasa = md.shrake_rupley(traj, probe_radius=0.14, mode='atom')[0]

    backbone_sasa = []
    sidechain_sasa = []
    seq = extract_sequence(pdb)

    for residue in traj.topology.residues:
        bb_indices = [atom.index for atom in residue.atoms if not atom.is_sidechain]
        sc_indices = [atom.index for atom in residue.atoms if atom.is_sidechain]

        bb_sasa = float(sum(sasa[i] for i in bb_indices)) * 100 #Correct from square nm to square angstroms
        sc_sasa = float(sum(sasa[i] for i in sc_indices)) * 100 #Correct from square nm to square angstroms

        backbone_sasa.append(bb_sasa)
        sidechain_sasa.append(sc_sasa)

    return seq, backbone_sasa, sidechain_sasa

def get_chain_info(pdb_path):
    """
    Retrieves amino acid sequences and backbone/sidechain SASA lists per chain,
    plus overall sequence and SASA for the entire protein.

    Parameters:
        pdb_path (str): Path to the PDB file.

    Returns:
        dict: keys = "Chain 1", "Chain 2", ..., "All"
            values = tuple (sequence:str or list, backbone_sasa:list, sidechain_sasa:list)
    """

    # Load structure
    traj = md.load_pdb(pdb_path)

    # Compute per-atom SASA once
    sasa_per_atom = md.shrake_rupley(traj, probe_radius=0.14, mode='atom')[0]  # 1D array

    # Extract sequences per chain from your helper (assumed returns list of strings)
    sequences = extract_sequences_by_chains(pdb_path)  # e.g. ["MKL...", "AGT..."]

    # Prepare result dictionary
    result = {}

    # MDTraj topology chains
    chains = list(traj.topology.chains)

    global_residue_index = 0  # Track residue position over all chains for "All"

    # Helper to sum SASA for atoms in residue by backbone/sidechain
    def residue_sasa(residue):
        bb_indices = [atom.index for atom in residue.atoms if not atom.is_sidechain]
        sc_indices = [atom.index for atom in residue.atoms if atom.is_sidechain]
        bb_s = float(sum(sasa_per_atom[i] for i in bb_indices)) * 100
        sc_s = float(sum(sasa_per_atom[i] for i in sc_indices)) * 100
        return bb_s, sc_s

    # Process each chain separately
    for i, chain in enumerate(chains, start=1):
        residues = list(chain.residues)
        backbone_sasa_list = []
        sidechain_sasa_list = []

        for res in residues:
            bb_s, sc_s = residue_sasa(res)
            backbone_sasa_list.append(bb_s)
            sidechain_sasa_list.append(sc_s)

        # Use the extracted sequence for the chain
        chain_seq = sequences[i - 1] if i-1 < len(sequences) else ""

        result[f"Chain {i}"] = (chain_seq, backbone_sasa_list, sidechain_sasa_list)

    # Compute overall (all residues in order)
    all_residues = list(traj.topology.residues)
    all_bb_sasa = []
    all_sc_sasa = []
    for res in all_residues:
        bb_s, sc_s = residue_sasa(res)
        all_bb_sasa.append(bb_s)
        all_sc_sasa.append(sc_s)

    # Get full sequence as one string using your helper
    full_seq = extract_sequence(pdb_path)

    result["All"] = (full_seq, all_bb_sasa, all_sc_sasa)

    return result

def sasa_to_rasa(seq, backbone_sasa, sidechain_sasa):
    """
    Converts backbone and sidechain SASA (Solvent Accessible Surface Area) values 
    to RASA (Relative Accessible Surface Area) values.

    Parameters:
        seq (str): The amino acid sequence.
        backbone_sasa (list): List of SASA values for the backbone.
        sidechain_sasa (list): List of SASA values for the sidechain.

    Returns:
        tuple: Two lists containing RASA values for the backbone and sidechain, respectively.
    """

    max_sasa = get_max_sasa_list()  # Retrieve the maximum SASA values for normalization

    backbone_max_sasa = max_sasa["backbone"]
    sidechain_max_sasa = max_sasa["sidechain"]  # Fixed: Previously referenced 'backbone' incorrectly

    # Compute RASA values by normalizing SASA with the corresponding maximum SASA values
    backbone_rasa = [backbone_sasa[i] / backbone_max_sasa[seq[i]] for i in range(len(seq))]
    sidechain_rasa = [sidechain_sasa[i] / sidechain_max_sasa[seq[i]] for i in range(len(seq))]

    return backbone_rasa, sidechain_rasa

def protein_unfolded_dG(pdb, osmolytes, custom_tfe=None, concentration=1.0, split_chains=False):
    """
    Computes the total free energy (dG) for the unfolded protein for one or multiple osmolytes.

    Parameters:
        pdb (str): Path to the PDB file.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        custom_tfe (dict): A dictionary with custom TFE values for osmolytes.
            WILL BE USED IN PLACE OF THE OSMOLYTE LIST
        concentration (float): The concentration in molar to scale the result. Default = 1.0.
        split_chains (bool): Whether or not the output should give a separate TFE value for each chain.

    Returns:
        dict: A dictionary where each key is an osmolyte (or chain, if split_chains = True) and the value is the total free energy.
    """
    # Ensure osmolytes is a list
    osmolytes = [osmolytes] if isinstance(osmolytes, str) else osmolytes

    # Choose sequence extraction method
    sequences = extract_sequences_by_chains(pdb) if split_chains else [extract_sequence(pdb)]
    full_sequence = extract_sequence(pdb) if split_chains else None

    results = {}
    for i, seq in enumerate(sequences, start=1):
        chain_results = {}

        chain_backbone_sasa, chain_sidechain_sasa = get_unfolded_sasa_from_sequence(seq)

        chain_backbone_rasa, chain_sidechain_rasa = sasa_to_rasa(seq, chain_backbone_sasa, chain_sidechain_sasa)

        for osmo in osmolytes:
            osmo_key = osmo.lower()

            chain_backbone_tfe, chain_sidechain_tfe = get_tfe(seq, osmo_key, custom_tfe)

            chain_backbone_result = concentration * np.sum(np.array([
                chain_backbone_rasa[i] * chain_backbone_tfe[i] for i in range(len(seq))]))
            chain_sidechain_result = concentration * np.sum(np.array([
                chain_sidechain_rasa[i] * chain_sidechain_tfe[i] for i in range(len(seq))]))

            chain_results[osmo_key] = chain_backbone_result + chain_sidechain_result

        if split_chains:
            results[f"Chain {i}"] = chain_results
        else:
            results.update(chain_results)

    if split_chains:
        chain_results = {}

        chain_backbone_sasa, chain_sidechain_sasa = get_unfolded_sasa_from_sequence(full_sequence)

        chain_backbone_rasa, chain_sidechain_rasa = sasa_to_rasa(full_sequence, chain_backbone_sasa, chain_sidechain_sasa)

        for osmo in osmolytes:
            osmo_key = osmo.lower()

            chain_backbone_tfe, chain_sidechain_tfe = get_tfe(full_sequence, osmo_key, custom_tfe)

            chain_backbone_result = concentration * np.sum(np.array([
                chain_backbone_rasa[i] * chain_backbone_tfe[i] for i in range(len(full_sequence))]))
            chain_sidechain_result = concentration * np.sum(np.array([
                chain_sidechain_rasa[i] * chain_sidechain_tfe[i] for i in range(len(full_sequence))]))

            chain_results[osmo_key] = chain_backbone_result + chain_sidechain_result
        results["All"] = chain_results

    return clean_dict(results)

def protein_folded_dG(pdb, osmolytes, custom_tfe=None, concentration=1.0, split_chains=False):
    """
    Computes the total free energy (dG) for the folded protein for one or multiple osmolytes.

    Parameters:
        pdb (str): Path to the PDB file.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        custom_tfe (dict, optional): A dictionary with custom TFE values for osmolytes. 
            WILL BE USED IN PLACE OF THE OSMOLYTE LIST
        concentration (float): The concentration in molar to scale the result. Default = 1.0.
        split_chains (bool): Whether or not the output should give a separate TFE value for each chain.

    Returns:
        dict: A dictionary where each key is an osmolyte (or chain, if split_chains = True), and the value is the total free energy.
    """
    # Ensure osmolytes is a list
    osmolytes = [osmolytes] if isinstance(osmolytes, str) else osmolytes

    # Choose sequence extraction method
    chains = get_chain_info(pdb) if split_chains else [get_pdb_info(pdb)]

    results = {}
    for i, chain in enumerate(chains, start=1):
        chain_results = {}
        protein = chains[chain] if split_chains else chain

        chain_seq, chain_backbone_sasa, chain_sidechain_sasa = protein
        chain_backbone_rasa, chain_sidechain_rasa = sasa_to_rasa(chain_seq, chain_backbone_sasa, chain_sidechain_sasa)

        for osmo in osmolytes:

            osmo_key = osmo.lower()

            chain_backbone_tfe, chain_sidechain_tfe = get_tfe(chain_seq, osmo_key, custom_tfe)

            chain_backbone_result = concentration * np.sum(np.array([
                chain_backbone_rasa[i] * chain_backbone_tfe[i] for i in range(len(chain_seq))]))
            chain_sidechain_result = concentration * np.sum(np.array([
                chain_sidechain_rasa[i] * chain_sidechain_tfe[i] for i in range(len(chain_seq))]))

            chain_results[osmo_key] = chain_backbone_result + chain_sidechain_result

        if split_chains and i != len(chains):
            results[f"Chain {i}"] = chain_results
        elif split_chains and i == len(chains):
            results["All"] = chain_results
        else:
            results.update(chain_results)

    return clean_dict(results)

def protein_ddG_folding(pdb, osmolytes, triplet=False, custom_tfe=None, concentration=1.0, split_chains = False):
    """
    Computes the change in free energy (dG) upon protein folding for one or multiple osmolytes.

    Parameters:
        pdb (str): Path to the PDB file.
        osmolytes (str or list): A single osmolyte or a list of osmolytes.
        triplet (bool): Whether to return the triplet (folded, unfolded, and their difference). Default = False.
        custom_tfe (dict, optional): A dictionary with custom TFE values for osmolytes.
            WILL BE USED IN PLACE OF THE OSMOLYTE LIST
        concentration (float): The concentration in molar to scale the result. Default = 1.0.
        split_chains (bool): Whether or not the output should give a separate TFE value for each chain.

    Returns:
        dict: A dictionary where each key is an osmolyte, and the value is either a tuple (folded_dG, unfolded_dG, dG_change) or the free energy difference.
    """
    # Ensure osmolytes is a list
    osmolytes = [osmolytes] if isinstance(osmolytes, str) else osmolytes

    # Compute folded and unfolded dG
    folded_dG = protein_folded_dG(pdb, osmolytes, custom_tfe, concentration, split_chains)
    unfolded_dG = protein_unfolded_dG(pdb, osmolytes, custom_tfe, concentration, split_chains)

    results = {}
    
    if split_chains:
        # Handle case where data is split by chains
        for chain in folded_dG.keys():  # Includes "All" as well
            results[chain] = {}
            for osmolyte in osmolytes:
                f_dG = folded_dG[chain].get(osmolyte, 0)  # Default to 0 if missing
                u_dG = unfolded_dG[chain].get(osmolyte, 0)
                ddG = f_dG - u_dG
                results[chain][osmolyte] = (f_dG, u_dG, ddG) if triplet else ddG
    else:
        # Handle case where data is not split by chains (just "All")
        for osmolyte in osmolytes:
            f_dG = folded_dG.get(osmolyte, 0)
            u_dG = unfolded_dG.get(osmolyte, 0)
            ddG = f_dG - u_dG
            results[osmolyte] = (f_dG, u_dG, ddG) if triplet else ddG
    
    return clean_dict(results)
