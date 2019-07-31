.. class:: center
.. image:: https://circleci.com/gh/usnistgov/jarvis.svg?style=shield
          :target: https://circleci.com/gh/usnistgov/jarvis
.. image:: https://travis-ci.org/usnistgov/jarvis.svg?branch=master
       :target: https://travis-ci.org/usnistgov/jarvis
.. image:: https://ci.appveyor.com/api/projects/status/d8na8vyfm7ulya9p/branch/master?svg=true
       :target: https://ci.appveyor.com/project/knc6/jarvis-63tl9
.. image:: https://api.codacy.com/project/badge/Grade/be8fa78b1c0a49c280415ce061163e77
       :target: https://www.codacy.com/app/knc6/jarvisutm_source=github.com&amp;utm_medium=referral&amp;utm_content=usnistgov/jarvis&amp;utm_campaign=Badge_Grade
.. image::  https://img.shields.io/pypi/dm/jarvis-tools.svg      
        :target: https://img.shields.io/pypi/dm/jarvis-tools.svg
.. image:: https://pepy.tech/badge/jarvis-tools
        :target: https://pepy.tech/badge/jarvis-tools
.. image:: https://codecov.io/gh/knc6/jarvis/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/knc6/jarvis

JARVIS
=====

Joint Automated Repository for Various Integrated Simulations (JARVIS) is an integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. The jarvis-tools package can be used for high-throughput computation, data-analysis, and training machine-learning models. Some of the packages used in the jarvis-tools package are shown below. JARVIS-official website: https://jarvis.nist.gov

.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/tools/jarvis-git.JPG
        :target: https://jarvis.nist.gov/
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/statistics.JPG
        :target: https://jarvis.nist.gov/
Installing JARVIS
-----------------
- We recommend installing miniconda environment from https://conda.io/miniconda.html ::

      bash Miniconda3-latest-Linux-x86_64.sh (for linux)
      bash Miniconda3-latest-MacOSX-x86_64.sh (for Mac)
      Download 32/64 bit python 3.6 miniconda exe and install (for windows)
      Now, let's make a conda environment just for JARVIS::
      conda create --name my_jarvis python=3.6
      source activate my_jarvis

- Git clone install (Recommended)::

      pip install numpy scipy matplotlib
      git clone https://github.com/usnistgov/jarvis.git
      cd jarvis
      python setup.py install


- Alternative pip install::

      pip install numpy scipy matplotlib
      pip install jarvis-tools

- Alternative nix install::
  Nix allows a robust and reproducible package for Linux. To generate a Nix environment for using JARVIS, follow the `Nix instructions`_.

.. _`Nix instructions`: ./nix/README.md

Jupyter notebooks
-----------------
- Python for beginners::
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/novice.JPG
        :target: https://colab.research.google.com/github/knc6/jarvis/blob/master/jarvis/colab/python_novice_notebook.ipynb
- JARVIS-DFT data analysis::
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/jdft.JPG
        :target: https://colab.research.google.com/github/knc6/jarvis/blob/master/jarvis/colab/jarvis_dft_explore_notebook.ipynb
- JARVIS-ML training::
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/jml_train.JPG
        :target: https://colab.research.google.com/github/knc6/jarvis/blob/master/jarvis/colab/jarvis_ml_quick_train_notebook.ipynb
- Comparing ML algorithms::
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/compareml.JPG
        :target: https://colab.research.google.com/github/knc6/jarvis/blob/master/jarvis/colab/compare_ml_algorithms_notebook.ipynb
- JARVIS-FF data-analysis::
.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/jff.JPG
        :target: https://colab.research.google.com/github/knc6/jarvis/blob/master/jarvis/colab/jarvis_ff_explore_notebook.ipynb
- See more in the plot-gallery below


References
-----------------
- JARVIS-FF::
      1) Evaluation and comparison of classical interatomic potentials through a user-friendly interactive web-interface, Nature: Sci Data. 4, 160125 (2017).https://www.nature.com/articles/sdata2016125
      2) High-throughput assessment of vacancy formation and surface energies of materials using classical force-fields, J. Phys. Cond. Matt. 30, 395901(2018).http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta
- JARVIS-DFT::
      3) High-throughput Identification and Characterization of Two-dimensional Materials using Density functional theory, Scientific Reports 7, 5179 (2017).https://www.nature.com/articles/s41598-017-05402-0
      4) Computational Screening of High-performance Optoelectronic Materials using OptB88vdW and TBmBJ Formalisms, Scientific Data 5, 180082 (2018).https://www.nature.com/articles/sdata201882
      5) Elastic properties of bulk and low-dimensional materials using van der Waals density functional, Phys. Rev. B, 98, 014107 (2018).https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.014107
      6) Convergence and machine learning predictions of Monkhorst-Pack k-points and plane-wave cut-off in high-throughput DFT calculations, Comp. Mat. Sci. 161, 300 (2019).https://www.sciencedirect.com/science/article/pii/S0927025619300813?via%3Dihub
      7) High-throughput Discovery of Topologically Non-trivial Materials using Spin-orbit Spillage, Nature: Sci. Rep. 9, 8534,(2019),  https://www.nature.com/articles/s41598-019-45028-y
      8) Accelerated Discovery of Efficient Solar-cell Materials using Quantum and Machine-learning Methods, Chem. Mater., https://pubs.acs.org/doi/10.1021/acs.chemmater.9b02166
      9) Data-driven Discovery of 3D and 2D Thermoelectric Materials , https://arxiv.org/abs/1903.06651.
- JARVIS-ML::
      10) Machine learning with force-field inspired descriptors for materials: fast screening and mapping energy landscape, Phys. Rev. Mat., 2, 083801 (2018).,https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
      11) Materials science in the artificial intelligence age: high-throughput library generation, machine learning, and a pathway from correlations to the underpinning physics, MRS Comm., 1-18 https://doi.org/10.1557/mrc.2019.95



Pypi, Readthedocs and Slideshare links
-----------------
      https://pypi.org/project/jarvis-tools
      
      https://jarvis-tools.readthedocs.io/en/latest/
      
      https://www.slideshare.net/KAMALCHOUDHARY4

Running the examples
-----------------
- For running high-throughput calculations, set HPC/system related information in env_variables
- Run py.test in tests folder to ensure basic setup
- LAMMPS example::
      An example calculation for Aluminum is given in the lammps folder for running EAM calculation (https://github.com/usnistgov/jarvis/blob/master/jarvis/lammps/examples/basic_input_output.py). Untar the example folder using tar -xvzf Al03.eam.alloy_nist.tgz . Change the 'parameters' variable and run jlammps.py.
- VASP example::
      Similarly, an example calculation for Silicon is given in vasp folder (https://github.com/usnistgov/jarvis/blob/master/jarvis/vasp/examples/runstruct_pyvasp.py). The input is a POSCAR file, which is already provided. executable paths, pseudopotential directory path and Special_POTCAR.yaml path needs to be adjusted in joptb88vdw.py top section. The master.py can be submitted to the queuing system with qsub sub.sh. 
- ML example::
      We trained machine learning models using JARVIS-DFT data on bandgaps, formation energies and elastic modulus and other properties. We used both chemical and structural descriptors during GradientBoostingRegression training. Example of getting 1557 descriptors for a system is given at: https://github.com/usnistgov/jarvis/blob/master/jarvis/sklearn/examples/desc_example.py
- Access to JARVIS database::
       Our database is freely available at https://www.ctcms.nist.gov/~knc6/JVASP.html, https://www.ctcms.nist.gov/jarvisml/, https://www.ctcms.nist.gov/~knc6/periodic.html, and https://www.ctcms.nist.gov/~knc6/JLAMMPS.html for JARVIS-DFT, JARVIS-ML and JARVIS-FF. 
       We can also load the dataset using python scripts similar to https://github.com/knc6/jarvis/blob/master/jarvis/db/static/explore_db.py .
- Uploading your data using JARVIS-API::
       In addition to downloading/browsing through the JARVIS-database, one can also upload their data and query using JARVIS-API. Follow the instructions in https://github.com/usnistgov/jarvis/blob/master/jarvis/db/mdcs/mdcs_api.py

Founders
-----------------
Kamal Choudhary, Francesca Tavazza (NIST)

Contributors
-----------------
Daniel Wheeler, Faical Yannick Congo, Kevin Garrity, Brian DeCost, Adam Biacchi,
Lucas Hale, Andrew Reid, Marcus Newrock (NIST)


Plot-gallery with additional jupyter notebooks
-----------------
.. class:: center
.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/RDF.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/RDF%2CPRDF%2CADF%2CDDF.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/ADF-a.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/RDF%2CPRDF%2CADF%2CDDF.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/ADF-b.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/RDF%2CPRDF%2CADF%2CDDF.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/DDF.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/RDF%2CPRDF%2CADF%2CDDF.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/bandstr.jpg
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/band_structure.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Dos.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/band_structure.ipynb

    
.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Wulff.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Wulff.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/BoltzTrap.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Boltztrap.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/kp_converg.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Convergence.ipynb

.. image:: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/en_converg.png
:Notebook: https://github.com/usnistgov/jarvis/blob/master/jarvis/db/static/Convergence.ipynb
