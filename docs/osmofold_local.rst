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
A list of maximum SASA values for each amino acid in the following order: `AFLIVPMGSTQNDEHKRCWY`.

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
The returned values represent the gTFE of the R-group ONLY. To get the combined gTFE of the R-group and backbone, append `"Back"` to the osmolyte name.  

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
- Trehalose (Hong *et al.*)

---

extract_sequences()
------------------

Extracts the amino acid sequence from a given PDB file as one-letter codes.

**Arguments:**

- **`pdb_file`**: A path to the PDB file of interest (relative to the working directory).  

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

- **`pdb_file`**: A path to the PDB file of interest (relative to the working directory).  

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

- **`seq`**: The amino acid sequence for which you want to compute TFE values.  

      Example: `"ACD"`

- **`osmo`**: The osmolyte you wish to compute with the given sequence.  

      Example: `"trehalose"`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom gTFE values, one for each of the 20 amino acids. Useful for testing osmolytes that OsmoFold doesn't currently support.

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

**Returns:**  
A list containing the gTFEs for a given sequence, with indices matching the amino acid sequence.

   Example: `[52.1, -31.2, 79.9]`

get_pdb_info()
--------

Returns the sequence and SASA for a given input PDB.

**Arguments:**

- **`pdb`**: A string containing the filepath of the input PDB.
      Example: `"/path/to/pdb.pdb"`

**Returns:**  
A list with two elements. [0] is the sequence of the input protein(s) as a string, and [1] is their corresponding SASA values.

   Example: `["ACD", [52.1, -31.2, 79.9]]`

get_chain_info()
--------

Returns the sequence and SASA for a given input PDB, split into individual chains.

**Arguments:**

- **`pdb`**: A string containing the filepath of the input PDB.

      Example: `"/path/to/pdb.pdb"`

**Returns:**  
A dictionary containing a key for each chain in the input PDB. Each corresponding value is a list with two elements, where 
[0] is the sequence of the input protein(s) as a string, and [1] is their corresponding SASA values. Also contains an "All"
key whose corresponding value will be the same as the output of get_pdb_info().

   Example: `{"Chain 1": ["ACD", [52.1, -31.2, 79.9]],
            "Chain 2": ["FPW", [-111.2, 90.4, 51.7]],
            "All": ["ACDFPW, [52.1, -31.2, 79.9, -111.2, 90.4, 51.7]]}`

sasa_to_rasa()
--------

Converts Solvent Accessible Surface Area (SASA) values into Relative Accessible Surface Area (RASA) values, where 1 
represents a fully exposed residue and 0 represents a fully buried residue.

**Arguments:**

- **`seq`**: The amino acid sequence for which you want to compute RASA values.  

      Example: `"ACD"`

- **`sasa_list`**: A list of SASA values with indices corresponding to the input sequence.

      Example: `[87.0, 135.2, 99.1]`

**Returns:**  
A list of RASA values with indices corresponding to the input sequence.

   Example: `[0.75, 0.89, 0.81]`

protein_unfolded_dG()
----------------------

Computes the total free energy (ΔG) for the unfolded protein in the presence of one or multiple osmolytes.

**Arguments:**

- **`pdb`**: The filepath to the input PDB file.  

      Example: `"/path/to/pdb.pdb"`

- **`osmolytes`**: A single osmolyte or a list of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`  
      Example: `["trehalose", "sucrose"]`

- **`backbone`**: OPTIONAL. Whether to include contributions from the protein backbone. Default is `True`.

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. The osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

- **`split_chains`**: OPTIONAL. Whether to compute ΔG separately for each protein chain. If `True`, the output will contain separate values for each chain. Default is `False`.

**Returns:**  
A dictionary where each key is an osmolyte (or a chain identifier if `split_chains=True`), and the corresponding value is the computed total free energy.

   Example (single-chain output):  

   `{"trehalose": -75.3, "sucrose": -42.1}`

   Example (multi-chain output with split_chains=True):

      `{
    "Chain 1": {"trehalose": -32.5, "sucrose": -18.4},
    "Chain 2": {"trehalose": -42.8, "sucrose": -23.7},
    "All": {"trehalose": -75.3, "sucrose": -42.1}
      }`

*If any of the functions fail to work as described, please submit a GitHub issue or contact Vincent (`vnichol2@uwyo.edu`).*
