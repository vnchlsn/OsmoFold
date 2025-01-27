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

- **`seq`**: The amino acid sequence you wish to run.  

      Example: `"ACD"`

- **`osmo`**: The osmolyte you wish to compute with the given sequence.  

      Example: `"trehalose"`

- **`custom_tfe`**: OPTIONAL. A dictionary of custom gTFE values, one for each of the 20 amino acid. Useful for testing osmolyte that OsmoFold doesn't currently support.

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

**Returns:**  
A list containing the gTFEs for a given sequence, with indicies matching the amino acid sequence.

   Example: `[52.1, -31.2, 79.9, ...]`

*If any of the functions fail to work as described, please submit a GitHub issue or contact Vincent (`vnichol2@uwyo.edu`).*
