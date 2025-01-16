.. A simple guide for quick calculations in osmofold

Usage
=====

Quick Start
-----------

This guide will help you quickly get started with OsmoFold. Use the script below to calculate the ΔΔG of osmolyte-induced protein folding and unfolding (commonly referred to as the "m-value").

**Note:** Before proceeding, ensure that you have installed OsmoFold correctly by referring to the :ref:`installation`.

.. code-block:: python

   from osmofold.osmofold_local import protein_ddG_folding

   output = protein_ddG_folding(
       pdb="path/to/pdb.pdb", 
       osmolytes=["trehalose", "tmao", "betaine"], 
       concentration=1.0,
       triplet = False
   )

   print(output)

### Explanation of Parameters

- **`pdb`**: Path to the input PDB file of the protein.
- **`osmolytes`**: List of osmolytes (e.g., "trehalose", "tmao", "betaine") to simulate their effects on protein folding.
- **`concentration`**: Concentration of the osmolyte(s) in molar (M). Default is `1.0` M.
- **`triplet`**: If `True`, the output will be a tuple containing the folded ΔG, unfolded ΔG, and ΔΔG. If `False`, only the ΔΔG (a float) will be returned. Default is `True`.

### Output

The results from `protein_ddG_folding` are reported in units of **calories per mole of protein (cal/mol)**:

- **Negative values**: Indicate stabilization of the **folded state**.
- **Positive values**: Indicate stabilization of the **unfolded state**.

Example Output:

.. code-block:: text

   {'trehalose': -50.0, 'tmao': -35.2, 'betaine': 12.5}

In this example, `trehalose` and `tmao` stabilize the folded state, while `betaine` stabilizes the unfolded state.

For advanced usage and additional functionality, consult the full documentation or reach out to the OsmoFold team.