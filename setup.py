.. class:: center
.. image:: https://circleci.com/gh/usnistgov/jarvis.svg?style=shield
    :target: https://circleci.com/gh/usnistgov/jarvis
.. image:: https://readthedocs.org/projects/jarvis-tools/badge/?version=latest&style=plastic
    :target: https://jarvis-tools.readthedocs.io/en/latest/?badge=latest
.. image:: https://api.codacy.com/project/badge/Grade/be8fa78b1c0a49c280415ce061163e77
    :target: https://www.codacy.com/app/knc6/jarvisutm_source=github.com&amp;utm_medium=referral&amp;utm_content=usnistgov/jarvis&amp;utm_campaign=Badge_Grade
.. image:: https://travis-ci.org/usnistgov/jarvis.svg?branch=master
    :target: https://travis-ci.org/usnistgov/jarvis



JARVIS
=====

Joint Automated Repository for Various Integrated Simulations (JARVIS) is an integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. Official website: https://jarvis.nist.gov



Installing JARVIS
-----------------
- We recommend installing miniconda environment from https://conda.io/miniconda.html .

      bash Miniconda3-latest-Linux-x86_64.sh (for linux)

      bash Miniconda3-latest-MacOSX-x86_64.sh (for Mac)

      Download 32/64 bit python 3.6 miniconda exe and install (for windows)

      Now, let's make a conda environment just for JARVIS::

      conda create --name my_jarvis python=3.6

      The 'my_jarvis' environment can be activated using the command

      source activate my_jarvis

- Simple pip install::

      pip install numpy scipy
      pip install jarvis-tools

- Git clone install::

      pip install numpy scipy
      git clone https://github.com/usnistgov/jarvis.git
      cd jarvis
      python setup.py install


- Nix allows a robust and reproducible package for Linux. To generate a Nix environment for using JARVIS, follow the `Nix instructions`_.

.. _`Nix instructions`: ./nix/README.md

Jupyter notebooks
-----------------
- JARVIS-FF:

      https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/jarvis_ff-explore.ipynb

- JARVIS-DFT:

      https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/jarvis_dft-explore.ipynb

- JARVIS-ML:

      https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/jarvis_ml-train.ipynb


How to cite JARVIS
-----------------
- JARVIS-FF:

      https://www.nature.com/articles/sdata2016125

      http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta

- JARVIS-DFT:

      https://www.nature.com/articles/s41598-017-05402-0

      https://www.nature.com/articles/sdata201882

      https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.014107

 - JARVIS-ML:

      https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801


Pypi
-----------------
- jarvis-tools:

      https://pypi.org/project/jarvis-tools

Readthedocs page
-----------------
- jarvis-tools:

      https://jarvis-tools.readthedocs.io/en/latest/

Running the examples
-----------------
- Set HPC/system related information in env_variables
- Run py.test in tests folder to ensure basic setup
- LAMMPS example: An example calculation for Aluminum is given in the lammps folder for running EAM calculation (https://github.com/usnistgov/jarvis/blob/master/jarvis/lammps/examples/basic_input_output.py). Untar the example folder using tar -xvzf Al03.eam.alloy_nist.tgz . Change the 'parameters' variable in run_alloy.py to accommodate your lammps executable path, inelast.mod and dsiplace.mod path according to your system. Right now, crystal structures are obtained from materials-project, so MAPI_KEY needs to be adjusted also. Then, type python run_alloy.py which will launch jobs on your cluster. A json file is created as any lammps calculation get created. First, the geometry is optimized using LAMMPS. Then, symmetrically distinct vacancies are created and their energetics is determined. Reference elements are also downloaded from MP and their energies are calculated for measuring the chemical potemtial used in defect calculation. Afterthat, surfaces with miller index upto 3 is created and their energies are calculated. Perfect cell is used as reference to calculate surface energies. Finally, phonons are calculated using phonopy package. We use ASE to make interface between phonopy and LAMMPS.
- VASP example: Similarly, an example calculation for Silicon is given in vasp folder (https://github.com/usnistgov/jarvis/blob/master/jarvis/vasp/examples/runstruct_pyvasp.py). The input is a POSCAR file, which is already provided. executable paths, pseudopotential directory path and Special_POTCAR.yaml path needs to be adjusted in joptb88vdw.py top section. The master.py can be submitted to the queuing system with qsub sub.sh. Adjust path and your cluster specifics in sub.sh. The contents of master.py and joptb88vdw.py should be independent of PBS/SLURM or any other cluster. After successful submission of the sub.sh script a series of ENCUT- and KPOINTS- files are produced for plane wave convergence and kpoint convergence. After that, geometry is optimized with force-convergence criteria in MAIN-RELAX file. The crystal structure used in MAIN-RELAX calculation is then used for subsequent band-structure, dielectric function, MBJ calculation and elastic property calculations. Each complete calculation is notified with generation of respective json file.
- ML example: We trained machine learning models using JARVIS-DFT data on bandgaps, formation energies and elastic modulus and other properties. We used both chemical and structural descriptors during GradientBoostingRegression training. Example of getting 1557 descriptors for a system is given at: https://github.com/usnistgov/jarvis/blob/master/jarvis/sklearn/examples/desc_example.py
- Access to JARVIS database: Our database is freely available at https://www.ctcms.nist.gov/~knc6/JVASP.html and https://www.ctcms.nist.gov/~knc6/periodic.html, https://www.ctcms.nist.gov/~knc6/JLAMMPS.html for JARVIS-DFT and JARVIS-FF. JARVIS-FF and its development version has data for energetics, elastic constants, defect formation energies, surface energies and phonon properties. JARVIS-DFT has 3D and 2D properties of materials such as energetics, X-ray diffraction (for 3D), radial-distribution function, electronic density of states, electrostatic potential and band-structure, MBJ and semi-local dielectric function calculation, elastic properties and allied gamma-point unconverged phonons and thermoelectric properties. JARVIS-DFT can be downloaded in json format also now using the link provided and JARVIS-DFT webpage.

Founders
-----------------
Kamal Choudhary, Francesca Tavazza (NIST)

Contributors
-----------------
Daniel Wheeler, Faical Yannick Congo, Kevin Garrity, Brian DeCost, Adam Biacchi,
Lucas Hale, Andrew Reid, Marcus Newrock (NIST)

Link to presentation slides
-----------------
https://www.slideshare.net/KAMALCHOUDHARY4
