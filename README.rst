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

NIST-JARVIS
=======================================

Joint Automated Repository for Various Integrated Simulations (JARVIS) is an integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. The jarvis-tools package can be used for high-throughput computation, data-analysis, and training machine-learning models. Some of the packages used in the jarvis-tools package are shown below. JARVIS-official website: https://jarvis.nist.gov

.. image:: https://github.com/knc6/jarvis/blob/master/jarvis/colab/colab_figures/statistics.JPG
        :target: https://jarvis.nist.gov/

Installing jarvis-tools
----------------------------------------

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

References
-----------------

- JARVIS-FF:
      1) Evaluation and comparison of classical interatomic potentials through a user-friendly interactive web-interface, Nature: Sci Data. 4, 160125 (2017).https://www.nature.com/articles/sdata2016125
      2) High-throughput assessment of vacancy formation and surface energies of materials using classical force-fields, J. Phys. Cond. Matt. 30, 395901(2018).http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta
- JARVIS-DFT:
      3) High-throughput Identification and Characterization of Two-dimensional Materials using Density functional theory, Scientific Reports 7, 5179 (2017).https://www.nature.com/articles/s41598-017-05402-0
      4) Computational Screening of High-performance Optoelectronic Materials using OptB88vdW and TBmBJ Formalisms, Scientific Data 5, 180082 (2018).https://www.nature.com/articles/sdata201882
      5) Elastic properties of bulk and low-dimensional materials using van der Waals density functional, Phys. Rev. B, 98, 014107 (2018).https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.014107
      6) Convergence and machine learning predictions of Monkhorst-Pack k-points and plane-wave cut-off in high-throughput DFT calculations, Comp. Mat. Sci. 161, 300 (2019).https://www.sciencedirect.com/science/article/pii/S0927025619300813?via%3Dihub
      7) High-throughput Discovery of Topologically Non-trivial Materials using Spin-orbit Spillage, Nature: Sci. Rep. 9, 8534,(2019),  https://www.nature.com/articles/s41598-019-45028-y
      8) Accelerated Discovery of Efficient Solar-cell Materials using Quantum and Machine-learning Methods, Chem. Mater., https://pubs.acs.org/doi/10.1021/acs.chemmater.9b02166
      9) Data-driven Discovery of 3D and 2D Thermoelectric Materials , https://arxiv.org/abs/1903.06651.
- JARVIS-ML:
      10) Machine learning with force-field inspired descriptors for materials: fast screening and mapping energy landscape, Phys. Rev. Mat., 2, 083801 (2018).,https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
      11) Materials science in the artificial intelligence age: high-throughput library generation, machine learning, and a pathway from correlations to the underpinning physics, MRS Comm., 1-18 https://doi.org/10.1557/mrc.2019.95



External links
-----------------------------------------
      https://pypi.org/project/jarvis-tools
      
      https://jarvis-tools.readthedocs.io/en/latest/
      
      https://www.slideshare.net/KAMALCHOUDHARY4

      https://figshare.com/authors/Kamal_Choudhary/4445539




