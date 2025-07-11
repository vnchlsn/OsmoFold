.. A simple guide for quick calculations in osmofold

Usage
=====

This guide will help you quickly get started with OsmoFold.

Quick Start
-----------

Use the script below to calculate the ΔΔG of osmolyte-induced protein folding and unfolding (commonly referred to as the "m-value").

**Note:** Before proceeding, ensure that you have installed OsmoFold correctly by referring to the `installation guide <https://osmofold.readthedocs.io/en/latest/installation.html>`_.

.. code-block:: python

   from osmofold.osmofold_local import protein_ddG_folding

   output = protein_ddG_folding(
       pdb="path/to/pdb.pdb", 
       osmolytes=["trehalose", "tmao", "betaine"], 
       concentration=1.0,
       triplet = False
   )

   print(output)

Explanation of Parameters

- **`pdb`**: Path to the input PDB file of the protein. This should be in the .pdb format (not .cif) and contain monomeric proteins with no non-protein atoms.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`concentration`**: Concentration of the osmolyte(s) in molar (M). Default is `1.0` M.
- **`triplet`**: If `True`, the output for each osmolyte will be a tuple containing the folded ΔG, unfolded ΔG, and ΔΔG. If `False`, only the ΔΔG (a float) will be returned. Default is `True`.

Output

The results from `protein_ddG_folding` are reported in units of **calories per mole of protein (cal/mol)**:

- **Negative values**: Indicate stabilization of the **folded state**.
- **Positive values**: Indicate stabilization of the **unfolded state**.

Example Output:

.. code-block:: text

   {'trehalose': -50.0, 'tmao': -35.2, 'betaine': 12.5}

In this example, `trehalose` and `tmao` stabilize the folded state (with trehalose doing so more robustly), while `betaine` stabilizes the unfolded state.

Individual TFEs by Chain
------------------------

The following approach applies for pdbs which contain multiple protein chains, if you want an individual TFE for each chains

.. code-block:: python

   from osmofold.osmofold_local import protein_ddG_folding

   output = protein_ddG_folding(
       pdb="path/to/pdb.pdb", 
       osmolytes=["trehalose", "tmao", "betaine"], 
       concentration=1.0,
       triplet = False,
       split_chains = True
   )

   print(output)

Explanation of Parameters

- **`pdb`**: Path to the input PDB file of the protein. This should be in the .pdb format (not .cif) and contain monomeric proteins with no non-protein atoms.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`concentration`**: Concentration of the osmolyte(s) in molar (M). Default is `1.0` M.
- **`triplet`**: If `True`, the output for each osmolyte will be a tuple containing the folded ΔG, unfolded ΔG, and ΔΔG. If `False`, only the ΔΔG (a float) will be returned. Default is `True`.
- **`split_chains`**: Determined whether results should be given for each peptide chain individually, just a sum for all chains. Default is `False`.

Output

The results from `protein_ddG_folding` are reported in units of **calories per mole of protein (cal/mol)**:

- **Negative values**: Indicate stabilization of the **folded state**.
- **Positive values**: Indicate stabilization of the **unfolded state**.

Example Output:

.. code-block:: text

   {'Chain 1': {'trehalose': -20.0, 'tmao': 10, 'betaine': 2.5},
   'Chain 2': {'trehalose': -30.0, 'tmao': -45.2, 'betaine': 10},
   'All': {'trehalose': -50.0, 'tmao': -35.2, 'betaine': 12.5}}

In this example, `trehalose` stabilizes the folded state of both chains. `tmao` stabilize the folded state of chain 2 but not chain 1, and `betaine` stabilizes the unfolded state of both chains.

Using Custom Solvent Accessible Surface Area (SASA) Values
------------------------

The default OsmoFold methods use MDTraj to calculate Solvent Accessible Surface Area (SASA) values.

However, OsmoFold does offer users the option to use SASA values calculated via other means. As a result, this additional functionality does not require a .pdb file to execute.

Below is an example of calculating the ΔΔG of a protein using custom SASA values:

.. code-block:: python

   from osmofold.osmofold_lite import protein_ddG_folding_lite

   output = protein_ddG_folding_lite(
       seq="SEQWENCE", 
       backbone_sasa = [92, 42, 38, 51, 82, 63, 71, 48]
       sidechain_sasa = [130, 98, 73, 76, 110, 115, 102, 180]
       osmolytes=["trehalose", "tmao", "betaine"], 
       concentration=1.0,
       triplet = False
   )

   print(output)

Explanation of Parameters

- **`seq`**: The one letter amino acid sequence for the protein of interest.
- **`backbone_sasa`**: A list contianing the SASA values for each residue in the protein sequence. These values should encompass both the backbone ONLY, with units of square angstroms.
- **`sidechain_sasa`**: A list contianing the SASA values for each residue in the protein sequence. These values should encompass both the sidechain ONLY, with units of square angstroms.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`concentration`**: Concentration of the osmolyte(s) in molar (M). Default is `1.0` M.
- **`triplet`**: If `True`, the output for each osmolyte will be a tuple containing the folded ΔG, unfolded ΔG, and ΔΔG. If `False`, only the ΔΔG (a float) will be returned. Default is `True`.

Output

The results from `protein_ddG_folding` are reported in units of **calories per mole of protein (cal/mol)**:

- **Negative values**: Indicate stabilization of the **folded state**.
- **Positive values**: Indicate stabilization of the **unfolded state**.

Example Output:

.. code-block:: text

   {'trehalose': -50.0, 'tmao': -35.2, 'betaine': 12.5}

In this example, `trehalose` and `tmao` stabilize the folded state (with trehalose doing so more robustly), while `betaine` stabilizes the unfolded state.

Batch Processing
----------------

Use the script below to batch process pdbs stored in the same folder. Unlike the previous functionality, it uses multithreading to handle large datasets more efficiently.

.. code-block:: python

   from osmofold.parallel import batch_process_pdbs

   if __name__ == "__main__":
      batch_process_pdbs(
         folder="path/to/pdbs", 
         osmolytes=["trehalose", "tmao", "betaine"], 
         save_csv = True,
         num_workers = 8,
         concentration = 1.0
      )
   
*Note: This `if __name__ == "__main__"` block is required for correct execution of batch processing.*

Explanation of Parameters

- **`folder`**: Path to the folder containing the pdbs you wish to test. As before, each should be in the .pdb format (not .cif) and contain monomeric proteins with no non-protein atoms.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`save_csv`**: Dictates whether to save the results to the local directory as a csv. Default is `True`.
- **`num_workers`**: The number of CPU cores that will be dedicated to running predictions. A safe bet for modern computers is 8. Default is `1`.
- **`concentration`**: Concentration of the osmolyte(s) in molar (M). Default is `1.0` M.

Output

For each protein-osmolyte combination, a new row will be created in the output csv, corresponding to the following columns.

- **`PDB name`**: The name of the pdb being tested.
- **`protein length`**: The length of the protein being tested.
- **`osmolyte`**: The osmolyte being tested.
- **`dG_Unfolded`**: The effect of the osmolyte on the free energy protein's unfolded state (assumes maximum solvent accessible surface area).
- **`dG_Folded`**: The effect of the osmolyte on the free energy protein's folded state (calculates solvent accessible surface areas from the provided pdb).
- **`ddG_Folding`**: The difference between ΔG_Folded and ΔG_Unfolded. Negative values indicate stabilization of the folded state. Positive values indicate stabilization of the unfolded state.
- **`error`**: Any errors that were thrown in the computation of that protein-osmolyte combination. For a successful prediction, should be `N/A`.

The results for all ΔG values are reported in units of **calories per mole of protein (cal/mol)**

*For advanced usage and additional functionality, consult the full documentation or reach out to Vincent (`vnichol2@uwyo.edu`).*
