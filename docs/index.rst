.. OsmoFold documentation master file

Overview
========

Vincent Nicholson¹ and Thomas C. Boothby¹

¹ Department of Molecular Biology, University of Wyoming, Laramie WY

**OsmoFold** is a tool that uses Auton and Bolen's empirical transfer free energy measurements to calculate the effect of osmolytes on protein folding. 

If you haven't already done so, we recommend reading the following literature:

   * Application of the transfer model to understand how naturally occurring osmolytes affect protein stability
   * Predicting the energetics of osmolyte-induced protein folding/unfolding
   * Additive transfer free energies of the peptide backbone unit that are independent of the model compound and the choice of concentration scale
   * Its Preferential Interactions with Biopolymers Account for Diverse Observed Effects of Trehalose

Below are some recent papers, preprints, and reviews that cover this method or an adapted version of it. They might be useful for your understanding of the protocol:

   * Disordered proteins interact with the chemical environment to tune their protective function during drying
   * LEA_4 motifs function alone and in conjunction with synergistic cosolutes to protect a labile enzyme during desiccation
   * Osmolyte-IDP Interactions During Desiccation

In short, OsmoFold will quantitatively assess the impact of several osmolytes on a given protein conformational change. Please note:

   * OsmoFold does not guarantee that said conformational change is physiologically relevant.
   * OsmoFold requires the existence of a known, binary conformational change in the protein(s) of interest.
   * OsmoFold cannot rule out the presence of specific protein-small molecule interactions, such as direct binding.

In general, **OsmoFold should be used to support experimental data, not replace it**.

General Installation
--------------------

To get started with OsmoFold, follow these steps:

   1. Set up and activate a Conda environment (Python\>=3.10):
         
      .. code-block:: python   
         
         conda create -n osmofold python=3.12

         conda activate osmofold

   2. Install numpy, mdtraj, and soursop:
      
      .. code-block:: python
         
         pip install numpy
         
         pip install mdtraj
         
         pip install soursop

   3. Pull and install OsmoFold:

      .. code-block:: python

         git pull https://github.com/vnchlsn/OsmoFold
         
         cd OsmoFold
         
         pip install .

   4. [Optional] Perform unit tests:

      .. code-block:: python

         pip install pytest

         pytest ./test

Version History 
---------------
WIP

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

   installation
   usage