JARVIS 
=====

NIST's Joint Automated Repository for Various Integrated Simulations (JARVIS) is a integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. JARVIS heavily uses VASP, LAMMPS, 
pymatgen, ase and scikit-learn packages.



Installing JARVIS
-----------------
- First we recommend installing miniconda environment from https://conda.io/miniconda.html .
  
      bash Miniconda3-latest-Linux-x86_64.sh (for linux)
      bash Miniconda3-latest-MacOSX-x86_64.sh (for Mac)
      Download 32/64 bit python 2.7 miniconda exe and install
- Now, let's make a conda environment just for JARVIS::

      conda create --name my_jarvis python=2.7
- The 'my_jarvis' environment can be activated using the command::

       source activate my_jarvis
       
- Then, get jarvis repo using the command::

      git clone git@github.com:usnistgov/jarvis.git .
- Go to the jarvis directory and type::

      python setup.py install      
 

How to cite JARVIS 
-----------------
- https://www.nature.com/articles/sdata2016125
- https://www.nature.com/articles/s41598-017-05402-0

-----------------
Author:Kamal Choudhary
