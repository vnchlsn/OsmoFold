.. _osmofold_local-docs:

**osmofold_local.py**
=====================

A summary of the functions in `osmofold_local.py` (the core script of OsmoFold) and how they work.

get_max_sasa_list()
------------------

A function that returns a list of maximum SASA (Solvent Accessible Surface Area) values for each amino acid.

**Arguments:**  
None.

**Returns:**  
A dictionary containing the unfolded SASA of protein sidechains and backbones.

.. code-block:: python

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
---

amino_to_energy()
--------------------

For a given osmolyte and a given amino acid, this function returns the experimentally derived gTFEs (Gibbs Transfer Free Energies) determined by Auton and Bolen.

**Arguments:**

- **`amino`**: The one-letter code for the amino acid you wish to get a gTFE for.

      Example: `'A'`

- **`cosolute`**: The osmolyte of interest.
   
      Example: `'trehalose'`

**Returns:**  
An integer representing the gTFE for the provided amino acid-osmolyte combination.

**Notes:**  
The returned values represent the gTFE of the R-group ONLY. To get the gTFE of the backbone only, append `"Back"` to the osmolyte name.

   Example: `"tmaoBack"`, `"trehaloseBack"`, `"ureaBack"`

**Supported Osmolytes**  
Each osmolyte has both R-group and backbone values:

- TMAO (Auton and Bolen)
- Sarcosine (Auton and Bolen)
- Betaine (Auton and Bolen)
- Sorbitol (Auton and Bolen)
- Sucrose (Auton and Bolen)
- Urea (Auton and Bolen)
- Proline (Auton and Bolen)
- Glycerol (Auton and Bolen)
- Trehalose (Auton and Bolen)
- Trehalose (Hong *et al.* 2015)

---

extract_sequences()
------------------

Extracts the amino acid sequence from a given PDB file as one-letter codes.

**Arguments:**

- **`pdb_file`**: A string containing the path to the PDB file of interest (relative to the working directory).  

      Example: `"your/path/here.pdb"`

**Returns:**  
A string containing the one-letter code for the protein in the specified PDB file.  

   Example:  `"SEQWENCE"`

**Notes:**  
This function is only compatible with PDB files containing protein chains. All other atoms must be removed.

---

extract_sequence_by_chains()
---------------------------

Extracts the amino acid sequence from a given PDB file as one-letter codes, divided by protein chain.

**Arguments:**

- **`pdb_file`**: A string containing the path to the PDB file of interest (relative to the working directory).  

      Example: `"your/path/here.pdb"`

**Returns:**  
A string containing the one-letter code for the individual protein chains in the specified PDB file. 

   Example:  `["SEQWENCE", "AAMTEARRD"]`

**Notes:**  
This function is only compatible with PDB files containing protein chains. All other atoms must be removed.

---

three_to_one()
-------------

A simple function for converting three-letter amino acid code to one-letter.

**Arguments:**

- **`residue`**: A string containing the three-letter amino acid code you wish to convert.  

      Example: `"ALA"`

**Returns:**  
A character containing the corresponding one-letter amino acid code.

   Example: `'A'`

---

get_tfe()
--------

Returns gTFEs for an entire protein sequence and a given osmolyte.

**Arguments:**

- **`seq`**: A string containing the amino acid sequence for which you want to compute TFE values.  

      Example: `"ACD"`

- **`osmo`**: A string containing the osmolyte you wish to compute with the given sequence.  

      Example: `"trehalose"`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom gTFE values, one for each of the 20 amino acids. This should cover both the backbone AND side chain.
Useful for testing osmolytes that OsmoFold doesn't currently support.

Each value key pair should be made up of a character (amino acid) and a float (gTFE).

      Example: 
      
      .. code-block:: python

            {
                "backbone": { 'A': X, 'F': Y, 'L': Z, ... },
                "sidechain": { 'A': A, 'F': A, 'L': A, ... }
            }

**Returns:**  
 Tuple: Two lists of TFE values for each amino acid in the sequence (backbone, sidechain).
 The indices correspond to the input sequence.

   Example: `([22, 22, 22], [52.1, -31.2, 79.9])`

---

get_pdb_info()
--------

Returns the sequence and SASA for a given input PDB.

**Arguments:**

- **`pdb`**: A string containing the filepath of the input PDB.
      Example: `"/path/to/pdb.pdb"`

**Returns:**  
 Tuple: Contians the protein sequence, and lists of SASA values for each amino acid in the sequence (seq, backbone, sidechain).
 The indices correspond to the input sequence.

   Example: `("ACD", [62.1, 55.2, 21.7], [33.1, 24.1, 19.7])`

---

get_chain_info()
--------

Returns the sequence and SASA for a given input PDB, split into individual chains.

**Arguments:**

- **`pdb`**: A string containing the filepath of the input PDB.

      Example: `"/path/to/pdb.pdb"`

**Returns:**  
A dictionary containing a key for each chain in the input PDB. Each corresponding value is a tuple with three elements, where 
the first is the sequence of the input protein(s) as a string, the second is their corresponding SASA values for the backbone 
stored as floats, and the 3rd is the corresponding SASA values for each sidechain (also stored as floats).

Also contains an "All" key whose corresponding value will be the same as the output of get_pdb_info().

   Example: 
   
   .. code-block:: python

            {"Chain 1": ("ACD", [62.1, 55.2, 21.7], [33.1, 24.1, 19.7]), 
            "Chain 2": ("FPW", [15.7, 21.6, 33.3], [65.1, 54.1, 41.2]), 
            "All": ("ACDFPW", [62.1, 55.2, 21.7, 15.7, 21.6, 33.3], [33.1, 24.1, 19.7, 65.1, 54.1, 41.2])}

---

sasa_to_rasa()
--------

Converts Solvent Accessible Surface Area (SASA) values into Relative Accessible Surface Area (RASA) values, where 1 
represents a fully exposed residue and 0 represents a fully buried residue.

**Arguments:**

- **`seq`**: A string containing the amino acid sequence for which you want to compute RASA values.  

      Example: `"ACD"`

- **`backbone_sasa`**: A list of SASA values for the backbone with indices corresponding to the input sequence, stored as floats.

      Example: `[62.1, 55.2, 21.7]`

- **`sidechain_sasa`**: A list of SASA values for the sidechains with indices corresponding to the input sequence, stored as floats.

      Example: `[87.0, 135.2, 99.1]`

**Returns:**  
A tuple containing two lists of RASA values for the backbone and sidechain.

   Example: `([0.8, 0.63, 0.21], [0.75, 0.43, 0.92])`

---

protein_unfolded_dG()
----------------------

Computes the total free energy (ΔG) for the unfolded protein in the presence of one or multiple osmolytes.

**Arguments:**

- **`pdb`**: A string containing the filepath to the input PDB file.  

      Example: `"/path/to/pdb.pdb"`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`

      Example: `["trehalose", "sucrose"]`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

- **`split_chains`**: OPTIONAL. A boolean that indicates whether to compute ΔG separately for each protein chain. If `True`, the output will contain separate values for each chain. Default is `False`.

**Returns:**  
A dictionary where each key is an osmolyte (or a chain identifier if `split_chains=True`), and the corresponding value is the computed total free energy.

   Example (single-chain output):  

      .. code-block:: python

            {"trehalose": -75.3, "sucrose": -42.1}

   Example (multi-chain output with split_chains=True):

     .. code-block:: python

            {
            "Chain 1": {"trehalose": -32.5, "sucrose": -18.4},
            "Chain 2": {"trehalose": -42.8, "sucrose": -23.7},
            "All": {"trehalose": -75.3, "sucrose": -42.1}
            }

---

protein_folded_dG()
----------------------

Computes the total free energy (ΔG) for the folded protein in the presence of one or multiple osmolytes.

**Arguments:**

- **`pdb`**: A string containing the filepath to the input PDB file.  

      Example: `"/path/to/pdb.pdb"`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`

      Example: `["trehalose", "sucrose"]`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

- **`split_chains`**: OPTIONAL. A boolean that indicates whether to compute ΔG separately for each protein chain. If `True`, the output will contain separate values for each chain. Default is `False`.

**Returns:**  
A dictionary where each key is an osmolyte (or a chain identifier if `split_chains=True`), and the corresponding value is the computed total free energy.

   Example (single-chain output):  

      .. code-block:: python
      
            {"trehalose": -53.7, "sucrose": -28.4}

   Example (multi-chain output with split_chains=True):

      .. code-block:: python

            {
            "Chain 1": {"trehalose": -21.4, "sucrose": -10.9},
            "Chain 2": {"trehalose": -32.3, "sucrose": -17.5},
            "All": {"trehalose": -53.7, "sucrose": -28.4}
            }

---

protein_ddG_folding()
----------------------

Computes the change in free energy (ΔΔG) of a protein conformational change for one or multiple osmolytes.

**Arguments:**

- **`pdb`**: A string containing the filepath to the input PDB file.  

      Example: `"/path/to/pdb.pdb"`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`  

      Example: `["trehalose", "sucrose"]`

- **`triplet`**: OPTIONAL. A boolean that determines whether the function returns a triplet containing the folded ΔG, unfolded ΔG, and their difference (ΔΔG). If `False`, only the free energy difference (ΔΔG) is returned. Default is `False`.

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

- **`split_chains`**: OPTIONAL. A boolean that indicates whether to compute ΔG separately for each protein chain. If `True`, the output will contain separate values for each chain. Default is `False`.

**Returns:**  
A dictionary where each key is an osmolyte (or a chain identifier if `split_chains=True`), and the corresponding value is either:  
- A floating-point value representing the free energy difference (ΔΔG).  
- A tuple `(folded_dG, unfolded_dG, ΔΔG)` if `triplet=True`.  

   Example (single-chain output with `triplet=False`):  

      .. code-block:: python
      
            {"trehalose": -22.5, "sucrose": -13.7}

   Example (single-chain output with `triplet=True`):  

      .. code-block:: python

            {"trehalose": (-53.7, -31.2, -22.5), "sucrose": (-28.4, -14.7, -13.7)}

   Example (multi-chain output with `split_chains=True` and `triplet=False`):

      .. code-block:: python
            
            {
            "Chain 1": {"trehalose": -10.1, "sucrose": -5.8},
            "Chain 2": {"trehalose": -12.4, "sucrose": -7.9},
            "All": {"trehalose": -22.5, "sucrose": -13.7}
            }

   Example (multi-chain output with `split_chains=True` and `triplet=True`):

      .. code-block:: python

            {
            "Chain 1": {"trehalose": (-21.4, -11.3, -10.1), "sucrose": (-10.9, -5.1, -5.8)},
            "Chain 2": {"trehalose": (-32.3, -19.9, -12.4), "sucrose": (-17.5, -9.6, -7.9)},
            "All": {"trehalose": (-53.7, -31.2, -22.5), "sucrose": (-28.4, -14.7, -13.7)}
            }

*If any of the functions fail to work as described, please submit a GitHub issue or contact Vincent (`vnichol2@uwyo.edu`).*
