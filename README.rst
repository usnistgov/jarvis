.. class:: center
.. image:: https://circleci.com/gh/usnistgov/jarvis.svg?style=shield
        :target: https://circleci.com/gh/usnistgov/jarvis
.. image:: https://travis-ci.org/usnistgov/jarvis.svg?branch=master
        :target: https://travis-ci.org/usnistgov/jarvis
.. image:: https://ci.appveyor.com/api/projects/status/d8na8vyfm7ulya9p/branch/master?svg=true
        :target: https://ci.appveyor.com/project/knc6/jarvis-63tl9 
.. image:: https://github.com/usnistgov/jarvis/workflows/JARVIS-Tools%20github%20action/badge.svg
        :target: https://github.com/usnistgov/jarvis
.. image:: https://github.com/usnistgov/jarvis/workflows/JARVIS-Tools%20linting/badge.svg
        :target: https://github.com/usnistgov/jarvis
.. image:: https://readthedocs.org/projects/jarvis-tools/badge/?version=latest
       :target: https://jarvis-tools.readthedocs.io/en/latest/?badge=latest   
.. image:: https://img.shields.io/codecov/c/github/knc6/jarvis
        :target: https://codecov.io/gh/knc6/jarvis  
.. image::  https://img.shields.io/pypi/dm/jarvis-tools.svg      
        :target: https://img.shields.io/pypi/dm/jarvis-tools.svg 
.. image:: https://pepy.tech/badge/jarvis-tools
        :target: https://pepy.tech/badge/jarvis-tools  
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3903515.svg
        :target: https://doi.org/10.5281/zenodo.3903515  
.. image:: https://img.shields.io/github/v/tag/usnistgov/jarvis
        :target: https://github.com/usnistgov/jarvis
.. image:: https://app.codacy.com/project/badge/Grade/be8fa78b1c0a49c280415ce061163e77    
        :target: https://www.codacy.com/manual/knc6/jarvis?utm_source=github.com&amp
.. image:: https://img.shields.io/github/commit-activity/y/usnistgov/jarvis   
        :target: https://github.com/usnistgov/jarvis
.. image:: https://img.shields.io/github/repo-size/usnistgov/jarvis   
        :target: https://github.com/usnistgov/jarvis
.. image:: https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Ftwitter.com%2Fjarvisnist
        :target: https://twitter.com/jarvisnist
.. image:: https://img.shields.io/badge/Facebook-Follow-Blue.svg
        :target: https://www.facebook.com/jarvisnist/
.. image:: https://img.shields.io/badge/LinkedIn-Follow-Blue.svg
        :target: https://www.linkedin.com/company/jarvisnist

        
========================================================================================

JARVIS-Tools: an open-source software package for data-driven atomistic materials design
=========================================================================================


NIST-JARVIS (Joint Automated Repository for Various Integrated Simulations) is an integrated framework for computational science using density functional theory,
classical force-field/molecular dynamics and machine-learning. The jarvis-tools package consists of scripts used in generating and analyzing the dataset. The NIST-JARVIS official website is: https://jarvis.nist.gov . This project is a part of the Materials Genome Initiative (MGI) at NIST (https://mgi.nist.gov/).

For more details, checkout our latest article:  `JARVIS: An Integrated Infrastructure for Data-driven Materials Design <https://arxiv.org/abs/2007.01831>`__

.. image:: https://www.ctcms.nist.gov/~knc6/images/logo/jarvis-mission.png
   :target: https://jarvis.nist.gov/



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
  ``JARVIS-WTB``      Wannier Tight Binding Hamiltonian parameter dataset. Website: https://www.ctcms.nist.gov/jarviswtb .
  ``JARVIS-EFG``      Electric field gradient dataset. Dataset will be made available and the website will be available soon.
  ``Downloads``       Download raw metadat at: https://www.ctcms.nist.gov/~knc6/downloads.html
  ===============  =======================================================================


Installation
---------------

Please see `Installation instructions <https://github.com/usnistgov/jarvis/blob/master/Installation.rst>`__


Example Jupyter notebooks
-----------------------------

Please find several `Google Colab Notebooks <https://github.com/JARVIS-Materials-Design/jarvis-tools-notebooks>`__


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

Please see `Publications related to JARVIS-Tools <https://github.com/usnistgov/jarvis/blob/master/Publications.rst>`__

External links
-----------------------------------------

      https://figshare.com/authors/Kamal_Choudhary/4445539
         
      https://pypi.org/project/jarvis-tools
      
      https://jarvis-tools.readthedocs.io/en/latest/
      
      https://www.slideshare.net/KAMALCHOUDHARY4

   


Correspondence
--------------------
Please report bugs as Github issues (https://github.com/usnistgov/jarvis/issues) or email to kamal.choudhary@nist.gov.

Funding support
--------------------

NIST-MGI (https://www.nist.gov/mgi).

