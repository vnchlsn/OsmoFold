.. A simple guide for quick calculations in osmofold

Usage
============

Quick Start
--------------------

The below script can be used to calculate the ΔΔG of osmolyte-induced protein folding and unfolding (commonly refered to as "m-value").

Be sure to refer to the "installation" page for a guide on how to install osmofold.

   .. code-block:: python

      from osmofold.osmofold_local import protein_ddG_folding

      output = protein_ddG_folding(pdb = "path/to/pdb.pdb", osmolytes = ["trehalose", "tmao", "betaine"], concentration = 1.0)

      print(output)

Here, negative results indicate stabilization of the folded state, and positive results indicate stabilization of the unfolded state. 

All values are returned with the units calories per mol of protein (cal/mol). 

The "concentration" argument specifies the concentration of the osmolyte of interest in molar, with the default being 1M.