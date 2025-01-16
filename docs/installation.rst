.. A simple guide for installing OsmoFold in it's current form

Installation
============

General Installation
--------------------

To install OsmoFold, follow these steps:

1. **Set up and activate a Conda environment (Python >= 3.10):**

   .. code-block:: bash

      conda create -n osmofold python=3.12

      conda activate osmofold

   *Note: The use of an environment manager such as Conda is not strictly required, but is highly reccomended. 
   I can not guarentee successful installation on an unclean environment.*

2. **Install dependencies:**

   Use the following commands to install the required libraries:

   .. code-block:: bash

      pip install numpy
      pip install mdtraj
      pip install soursop

3. **Clone and install OsmoFold:**

   Download the OsmoFold repository and install it:

   .. code-block:: bash

      git clone https://github.com/vnchlsn/OsmoFold
      cd OsmoFold
      pip install .

4. **[Optional] Run unit tests:**

   To verify your installation, you can run the provided unit tests:

   .. code-block:: bash

      pip install pytest
      cd test
      pytest

   
**Note:** The above instructions have been tested on multiple Unix-based systems and Windows (using WSL2). 
If you encounter any issues, please email `vnichol2@uwyo.edu` or submit a GitHub issue.

Other Installation Methods
--------------------------

Currently, OsmoFold does **not** support installation via the following methods:

1. **Conda:**

   .. code-block:: bash

      conda install osmofold

2. **Poetry:**

   .. code-block:: bash

      poetry install osmofold

3. **Pipenv:**

   .. code-block:: bash

      pipenv install osmofold

Thank you for using OsmoFold! If you have any questions or feedback, feel free to reach out.