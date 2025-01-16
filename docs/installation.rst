.. A simple guide for installing OsmoFold in it's current form

Installation
============

General Installation
--------------------

To get started with OsmoFold, follow these steps:

   1. Set up and activate a Conda environment (Python\>=3.10):
         
      .. code-block:: bash
         
         conda create -n osmofold python=3.12

         conda activate osmofold

   2. Install numpy, mdtraj, and soursop:
      
      .. code-block:: bash
         
         pip install numpy
         
         pip install mdtraj
         
         pip install soursop

   3. Pull and install OsmoFold:

      .. code-block:: bash

         git clone https://github.com/vnchlsn/OsmoFold
         
         cd OsmoFold
         
         pip install .

   4. [Optional] Perform unit tests:

      .. code-block:: bash

         pip install pytest

         cd test

         pytest

The above instructions have been tested on multiple versions of unix and windows (using WSL2). If the above does NOT work please email vnichol2@uwyo.edu or submit a github issue.

Other Methods
-------------

OsmoFold does NOT currently support installation via:

    1. Conda

        .. code-block:: bash
            
         conda install osmofold

    2. Poetry

        .. code-block:: bash

         poetry install

    3. Pipenv

        .. code-block:: bash
        
         pipenv install


