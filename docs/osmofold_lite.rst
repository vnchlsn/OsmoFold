.. _osmofold_lite-docs:

**osmofold_lite.py**
=====================

A summary of the functions in `osmofold_lite.py` (designed to perform TFE calculations without SOURSOP) and how they work.

read_fasta()
--------------------

Reads a FASTA file and returns the sequences within as a list of strings.

**Arguments:**

- **`fasta_path`**: A string containing the path to the FASTA file of interest (relative to the working directory).

**Returns:**  
A list of strings containing the sequences within the provided FASTA.

---

protein_unfolded_dG_lite()
----------------------

Computes the total free energy (ΔG) for the unfolded protein in the presence of one or multiple osmolytes.

**Arguments:**

- **`seq`**: A string containing the input sequence, as one-letter amino acid code.

      Example: `"ACD"`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`

      Example: `["trehalose", "sucrose"]`

- **`backbone`**: OPTIONAL. A boolean that indicates whether to include contributions from the protein backbone. Default is `True`.

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

**Returns:**  
A dictionary where each key is an osmolyte and the corresponding value is the computed total free energy.

   Example:  

      .. code-block:: python

            {"trehalose": -75.3, "sucrose": -42.1}

---

protein_folded_dG_lite()
----------------------

Computes the total free energy (ΔG) for the folded protein in the presence of one or multiple osmolytes.

**Arguments:**

- **`seq`**: A string containing the input sequence, as one-letter amino acid code.

      Example: `"ACD"`

- **`sasa_list`**: A list of floating-point values containing SASA values with indecies correpsonding to the input sequence.

      Example: `[87.0, 135.2, 99.1]`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`

      Example: `["trehalose", "sucrose"]`

- **`backbone`**: OPTIONAL. A boolean that indicates whether to include contributions from the protein backbone. Default is `True`.

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

**Returns:**  
A dictionary where each key is an osmolyte and the corresponding value is the computed total free energy.

   Example:  

      .. code-block:: python

            {"trehalose": -52.1, "sucrose": -27.4}

---

protein_ddG_folding()
----------------------

Computes the change in free energy (ΔΔG) of a protein conformational change for one or multiple osmolytes.

**Arguments:**

- **`seq`**: A string containing the input sequence, as one-letter amino acid code. 

      Example: `"ACD"`

- **`sasa_list`**: A list of floating-point values containing SASA values with indecies correpsonding to the input sequence.

      Example: `[87.0, 135.2, 99.1]`

- **`osmolytes`**: A string containing a single osmolyte or a list of strings of osmolytes to compute ΔG values for.  

      Example: `"trehalose"`  

      Example: `["trehalose", "sucrose"]`

- **`backbone`**: OPTIONAL. A boolean that indicates whether to include contributions from the protein backbone. Default is `True`.

- **`triplet`**: OPTIONAL. A boolean that determines whether the function returns a triplet containing the folded ΔG, unfolded ΔG, and their difference (ΔΔG). If `False`, only the free energy difference (ΔΔG) is returned. Default is `False`.

- **`custom_tfe`**: OPTIONAL. A dictionary of custom transfer free energy (TFE) values for specific osmolytes.  

      Example: `{'A': 52.1, 'C': -31.2, 'D': 79.9, ...}`

- **`concentration`**: OPTIONAL. A floating-point value denoting the osmolyte concentration in molar (M) to scale the computed free energy. Default is `1.0`.

**Returns:**  
A dictionary where each key is an osmolyte and the corresponding value is either:  
- A floating-point value representing the free energy difference (ΔΔG).  
- A tuple `(folded_dG, unfolded_dG, ΔΔG)` if `triplet=True`.  

   Example (single-chain output with `triplet=False`):  

      .. code-block:: python
      
            {"trehalose": -22.5, "sucrose": -13.7}

   Example (single-chain output with `triplet=True`):  

      .. code-block:: python

            {"trehalose": (-53.7, -31.2, -22.5), "sucrose": (-28.4, -14.7, -13.7)}

*If any of the functions fail to work as described, please submit a GitHub issue or contact Vincent (`vnichol2@uwyo.edu`).*
