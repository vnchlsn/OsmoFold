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

Batch Processing
----------------

Use the script below to batch process pdbs stored in the same folder. Unlike the previous functionality, it uses multithreading to handle large datasets more efficiently.

.. code-block:: python

   from osmofold.parallel import batch_process_pdbs

   if __name__ == "__main__":
      batch_process_pdbs(
         folder="path/to/pdbs", 
         osmolytes=["trehalose", "tmao", "betaine"], 
         save_csv = True
         num_workers = 8
      )

Explanation of Parameters

- **`folder`**: Path to the folder containing the pdbs you with to test. As before, each should be in the .pdb format (not .cif) and contain monomeric proteins with no non-protein atoms.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`save_csv`**: Dictates whether to save the results to the local directory as a csv. Default is `True`.
- **`num_workers`**: The number of CPU cores that will be dedicated to running predictions. A safe bet for modern computers is 8. Default is `1`.

Output

For each protien-osmolyte combination, a new row will be created in the output csv, corresponding to the following columns.

- **`PDB name`**: The name of the pdb being tested.
- **`protein length`**: The length of the protein being tested.
- **`osmolyte`**: The osmolyte being tested.
- **`dG_Unfolded`**: The effect of the osmolyte on the free energy protein's unfolded state (assumes maximum solvent accessible surface area).
- **`dG_Folded`**: The effect of the osmolyte on the free energy protein's folded state (calculates solvent accessible surface areas from the provided pdb).
- **`ddG_Folding`**: The difference between dG_Folded and dG_Unfolded. Negative values indicate stabilization of the folded state. Positive values indicate stabilization of the unfolded state.
- **`error`**: Any errors that were thrown in the computation of that protein-osmolyte combination. For a successful prediction, should be `N/A`.

*For advanced usage and additional functionality, consult the full documentation or reach out to Vincent (`vnichol2@uwyo.edu`).*