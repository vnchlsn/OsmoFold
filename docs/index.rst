.. OsmoFold documentation master file

Overview
========

An Introduction to OsmoFold
---------------------------

Vincent Nicholson¹ and Thomas C. Boothby¹

¹ Department of Molecular Biology, University of Wyoming, Laramie WY

**OsmoFold** is a tool that uses Auton and Bolen's empirical transfer free energy measurements to calculate the effect of osmolytes on protein folding. 

If you haven't already done so, we recommend reading the following literature:

   * `Application of the transfer model to understand how naturally occurring osmolytes affect protein stability <https://pubmed.ncbi.nlm.nih.gov/17875431/>`_
   * `Predicting the energetics of osmolyte-induced protein folding/unfolding <https://pubmed.ncbi.nlm.nih.gov/16214887/>`_
   * `Additive transfer free energies of the peptide backbone unit that are independent of the model compound and the choice of concentration scale <https://pubmed.ncbi.nlm.nih.gov/14756570/>`_
   * `Its Preferential Interactions with Biopolymers Account for Diverse Observed Effects of Trehalose <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4572414/>`_

Below are some recent papers, preprints, and reviews that cover this method or an adapted version of it. They might be useful for your understanding of the protocol:

   * `Disordered proteins interact with the chemical environment to tune their protective function during drying <https://elifesciences.org/articles/97231>`_
   * `LEA_4 motifs function alone and in conjunction with synergistic cosolutes to protect a labile enzyme during desiccation <https://www.biorxiv.org/content/10.1101/2024.09.04.611296v1.full.pdf>`_
   * `Osmolyte-IDP Interactions During Desiccation <https://www.sciencedirect.com/science/article/pii/S1877117324001765?via%3Dihub>`_

In short, OsmoFold will quantitatively assess the impact of several osmolytes on a given protein conformational change. Please note:

   * OsmoFold does not guarantee that said conformational change is physiologically relevant.
   * OsmoFold requires the existence of a known, binary conformational change in the protein(s) of interest.
   * OsmoFold cannot rule out the presence of specific protein-small molecule interactions, such as direct binding.

In general, **OsmoFold should be used to support experimental data, not replace it**.

Version History
---------------

* 1/23/25 - OsmoFold v0.3.0

Added the ability to give individual results for multiple interacting chains within the same pdb (see `usage <https://osmofold.readthedocs.io/en/latest/usage.html>`_).

* 1/16/25 - OsmoFold v0.2.2

Public release of unit tests for all current functionality

* 1/8/25 - OsmoFold v0.2.1

Added support for ΔG calculations from sequence with user-provided SASA values

* 12/30/24 - OsmoFold v0.2.0

Added support for mulithreading / batch processing of monomeric proteins

* 06/18/24 - OsmoFold v0.1.0

First public release. Support for basic ΔG calculations of monomeric proteins (PDB only)

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

   installation
   usage
   osmofold_local.py