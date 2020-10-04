.. class:: center
.. image:: https://circleci.com/gh/usnistgov/jarvis.svg?style=shield
        :target: https://circleci.com/gh/usnistgov/jarvis
.. image:: https://img.shields.io/travis/usnistgov/jarvis/master.svg?label=Travis%20CI
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


Some important features
=======================================================================

- **Software workflow tasks**:  VASP, Quantum Espresso, Wien2k BoltzTrap, Wannier90, LAMMPS, Scikit-learn, TensorFlow, LightGBM.

- **HPC clusters**: PBS and SLURM.

- **Examples**: Notebooks and test scripts to explain the package.

- **Available datasets**: `Summary of several datasets <https://github.com/usnistgov/jarvis/blob/master/DatasetSummary.rst>`__ .


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
>>> from jarvis.io.vasp.inputs import Poscar
>>> for i in dft_3d:
...     atoms = Atoms.from_dict(i['atoms'])
...     poscar = Poscar(atoms)
...     jid = i['jid']
...     filename = 'POSCAR-'+jid+'.vasp'
...     poscar.write_file(filename)
>>> dft_2d = data(dataset='dft_2d')
>>> print (len(dft_2d))
1070
>>> for i in dft_2d:
...     atoms = Atoms.from_dict(i['atoms'])
...     poscar = Poscar(atoms)
...     jid = i['jid']
...     filename = 'POSCAR-'+jid+'.vasp'
...     poscar.write_file(filename)

References
-----------------

Please see `Publications related to JARVIS-Tools <https://github.com/usnistgov/jarvis/blob/master/Publications.rst>`__

Documentation
-----------------------------------------
      https://jarvis-materials-design.github.io/dbdocs/



Correspondence
--------------------
Please report bugs as Github issues (https://github.com/usnistgov/jarvis/issues) or email to kamal.choudhary@nist.gov.

Funding support
--------------------

NIST-MGI (https://www.nist.gov/mgi).

