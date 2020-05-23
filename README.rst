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

.. image:: https://codecov.io/gh/usnistgov/jarvis/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/usnistgov/jarvis      
.. image:: https://www.ctcms.nist.gov/~knc6/jlogo.png
        :target: https://jarvis.nist.gov/

        
========================================================================================

jarvis-tools: an open-source software package for data-driven atomistic materials design
=========================================================================================




NIST-JARVIS (Joint Automated Repository for Various Integrated Simulations) is an integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. The jarvis-tools package consists of scripts used in generating and analyzing the dataset. The NIST-JARVIS official website is: https://jarvis.nist.gov . This project is a part of the Materials Genome Initiative (MGI) at NIST (https://mgi.nist.gov/).

* A summary of the projects

  ===============  =======================================================================
  Projects          Brief description
  ===============  =======================================================================
  ``JARVIS-DFT``      Density functional theory calculation database for ~40000 3D and ~1000 2D materials. Some of the material-properties include: Heat of formation, Crystal-structural data using OptB88vdW, PBE, LDA functionals, Bandgaps using semi-local, meta-GGA, HSE06 and other beyond DFT methods, Electron and phonon-bandstructures, Elastic, Piezoelectric, Thermoelectric, Dielectric tensors, Exfoliation energies for low-diemnsional materials, Frequency dependent dielectric function, Absorption coefficients, Work-function for 2D materials, Infrared and Raman intensities, Electric field gradient, Magnetic moment, Solar-cell efficiencies, Scanning Tunneling Microscopy (STM) images, Topological spin-orbit spillage, converged k-point and plane wave cut-offs, Wannier-tight binding Hamiltonian parameters and more. The website for JARVIS-DFT: https://www.ctcms.nist.gov/~knc6/JVASP.html
  ``JARVIS-FF``       Classical molecular dynamics calculation database for ~2000 3D materials with interatomic potential/force-fields. Some of the properties included in JARVIS-FF are energetics, elastic constants, surface energies, defect formations energies and phonon frequencies of materials. The website for JARVIS-FF: https://www.ctcms.nist.gov/~knc6/periodic.html
  ``JARVIS-ML``       Machine learning prediction tools trained on the JARVIS-DFT data. Some of the ML-prediction models are for  Heat of formation, GGA/METAGGA bandgaps, Refractive indices, Bulk and shear modulus, Magnetic moment, Thermoelectric, Piezoelectric and Dielectric properties properties, Exfoliation energies, Solar-cell efficiency, and STM image classification. The website for JARVIS-ML: https://www.ctcms.nist.gov/jarvisml/
  ``JARVIS-Het.``     Heterostructure design tools for 2D materials in the JARVIS-DFT database. Some of the properties available are: work function, Band-alignment, and Heterostructure classification. JARVIS-Heterostructure website: https://www.ctcms.nist.gov/jarvish/
  ``JARVIS-PV``       Solar-cell/Photovoltaic cell design tools. Dataset is made available and the website will be available soon.
  ``JARVIS-STM``      Scanning-tunneling microscopy images for 2D materials. Dataset is made available and the website will be available soon.
  ``JARVIS-WTB``      Wannier Tight Binding Hamiltonian parameter dataset. Dataset will be made available and the website will be available soon.
  ``JARVIS-EFG``      Electric field gradient dataset. Dataset will be made available and the website will be available soon.
  ===============  =======================================================================



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

Example Jupyter notebooks
-----------------------------
Look into the notebooks folder

Example function
-----------------
>>> from jarvis.core.atoms import Atoms
>>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
>>> coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
>>> elements = ["Si", "Si"]
>>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
>>> density = round(Si.density,2)
>>> print (density)
2.33
>>>
>>> from jarvis.db.figshare import data
>>> dft_3d = data(dataset='dft_3d')
>>> print (len(dft_3d))
36099

References
-----------------
- JARVIS-FF:
      1) Evaluation and comparison of classical interatomic potentials through a user-friendly interactive web-interface, Nature: Sci Data. 4, 160125 (2017). https://www.nature.com/articles/sdata2016125
      2) High-throughput assessment of vacancy formation and surface energies of materials using classical force-fields, J. Phys. Cond. Matt. 30, 395901(2018). http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta

- JARVIS-DFT:
      3) High-throughput Identification and Characterization of Two-dimensional Materials using Density functional theory, Scientific Reports 7, 5179 (2017). https://www.nature.com/articles/s41598-017-05402-0
      4) Computational Screening of High-performance Optoelectronic Materials using OptB88vdW and TBmBJ Formalisms, Scientific Data 5, 180082 (2018). https://www.nature.com/articles/sdata201882
      5) Elastic properties of bulk and low-dimensional materials using van der Waals density functional, Phys. Rev. B, 98, 014107 (2018). https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.014107
      6) High-throughput Discovery of Topologically Non-trivial Materials using Spin-orbit Spillage, Nature: Sci. Rep. 9, 8534,(2019), https://www.nature.com/articles/s41598-019-45028-y
      7) Computational Search for Magnetic and Non-magnetic 2D Topological Materials using Unified Spin-orbit Spillage Screening, npj Comp. Mat., 6, 49 (2020). https://www.nature.com/articles/s41524-020-0319-4 .
 

- JARVIS-ML:
      8) Machine learning with force-field inspired descriptors for materials: fast screening and mapping energy landscape, Phys. Rev. Mat., 2, 083801 (2018). https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
      9) Convergence and machine learning predictions of Monkhorst-Pack k-points and plane-wave cut-off in high-throughput DFT calculations, Comp. Mat. Sci. 161, 300 (2019). https://www.sciencedirect.com/science/article/pii/S0927025619300813?via%3Dihub
      10) Materials science in the artificial intelligence age: high-throughput library generation, machine learning, and a pathway from correlations to the underpinning physics, MRS Comm., 1-18, 2019. https://doi.org/10.1557/mrc.2019.95
      11) Enhancing materials property prediction by leveraging computational and experimental data using deep transfer learning, Nature Comm., 10, 1, (2019). https://www.nature.com/articles/s41467-019-13297-w
      12) Accelerated Discovery of Efficient Solar-cell Materials using Quantum and Machine-learning Methods, Chem. Mater., https://pubs.acs.org/doi/10.1021/acs.chemmater.9b02166
      13) High-throughput Density Functional Perturbation Theory and Machine Learning Predictions of Infrared, Piezoelectric and Dielectric Responses, https://arxiv.org/abs/1910.01183.
      14) Data-driven Discovery of 3D and 2D Thermoelectric Materials , https://arxiv.org/abs/1903.06651.

External links
-----------------------------------------
      https://pypi.org/project/jarvis-tools
      
      https://jarvis-tools.readthedocs.io/en/latest/
      
      https://www.slideshare.net/KAMALCHOUDHARY4

      https://figshare.com/authors/Kamal_Choudhary/4445539


Correspondence
--------------------
Please report bugs as Github issues (https://github.com/usnistgov/jarvis/issues) or email to kamal.choudhary@nist.gov.

Funding support
--------------------

NIST-MGI (https://www.nist.gov/mgi).

