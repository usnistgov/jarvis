Installing jarvis-tools
----------------------------------------

- We recommend installing miniconda environment from https://conda.io/miniconda.html ::

      bash Miniconda3-latest-Linux-x86_64.sh (for linux)
      bash Miniconda3-latest-MacOSX-x86_64.sh (for Mac)
      Download 32/64 bit python 3.6 miniconda exe and install (for windows)
      Now, let's make a conda environment just for JARVIS::
      conda create --name my_jarvis python=3.8
      source activate my_jarvis

- Git clone install (Recommended)::

      pip install numpy scipy matplotlib
      git clone https://github.com/usnistgov/jarvis.git
      cd jarvis
      python setup.py install


- Alternative pip install::

      pip install numpy scipy matplotlib
      pip install jarvis-tools

- Alternative conda install::
 
      conda install -c conda-forge jarvis-tools
