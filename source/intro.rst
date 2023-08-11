JARVIS-Tools
=========================================================================================

.. image:: https://badge.fury.io/py/jarvis-tools.svg
        :target: https://pypi.org/project/jarvis-tools/ 
.. image:: https://anaconda.org/conda-forge/jarvis-tools/badges/version.svg
        :target: https://anaconda.org/conda-forge/jarvis-tools   
.. image:: https://img.shields.io/github/v/tag/usnistgov/jarvis
        :target: https://github.com/usnistgov/jarvis
.. image:: https://img.shields.io/travis/usnistgov/jarvis/master.svg?label=Travis%20CI
        :target: https://travis-ci.org/usnistgov/jarvis
.. image:: https://ci.appveyor.com/api/projects/status/d8na8vyfm7ulya9p/branch/master?svg=true
        :target: https://ci.appveyor.com/project/knc6/jarvis-63tl9 
.. image:: https://github.com/usnistgov/jarvis/workflows/JARVIS-Tools%20github%20action/badge.svg
        :target: https://github.com/usnistgov/jarvis
.. image:: https://github.com/usnistgov/jarvis/workflows/JARVIS-Tools%20linting/badge.svg
        :target: https://github.com/usnistgov/jarvis  
.. image:: https://img.shields.io/codecov/c/github/knc6/jarvis
        :target: https://codecov.io/gh/knc6/jarvis  
.. image::  https://img.shields.io/pypi/dm/jarvis-tools.svg      
        :target: https://img.shields.io/pypi/dm/jarvis-tools.svg 
.. image:: https://pepy.tech/badge/jarvis-tools
        :target: https://pepy.tech/badge/jarvis-tools  
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3903515.svg
        :target: https://doi.org/10.5281/zenodo.3903515  
.. image:: https://app.codacy.com/project/badge/Grade/be8fa78b1c0a49c280415ce061163e77    
        :target: https://www.codacy.com/manual/knc6/jarvis?utm_source=github.com&amp
.. image:: https://img.shields.io/github/commit-activity/y/usnistgov/jarvis   
        :target: https://github.com/usnistgov/jarvis
.. image:: https://img.shields.io/github/repo-size/usnistgov/jarvis   
        :target: https://github.com/usnistgov/jarvis
.. image:: https://img.shields.io/badge/JARVIS-Figshare-Green.svg  
        :target: https://figshare.com/authors/Kamal_Choudhary/4445539
.. image:: https://img.shields.io/badge/JARVIS-DBDocs-Green.svg  
        :target: https://jarvis-materials-design.github.io/dbdocs   
.. image:: https://img.shields.io/badge/JARVIS-ToolsDocs-Green.svg  
        :target: https://jarvis-tools.readthedocs.io/en/latest/index.html 
.. image:: https://readthedocs.org/projects/jarvis-tools/badge/?version=master
       :target: https://jarvis-tools.readthedocs.io/en/master/?badge=latest  
.. image:: https://colab.research.google.com/assets/colab-badge.svg  
       :target: https://github.com/JARVIS-Materials-Design/jarvis-tools-notebooks



About
---------------------

The JARVIS-Tools is an open-access software package for atomistic data-driven materials desgin. JARVIS-Tools can be used for:
a) setting up calculations, b) analysis and informatics, c) plotting, d) database development, e) machine-learning, and f) web-page development.

JARVIS-Tools empowers NIST-JARVIS (Joint Automated Repository for Various Integrated Simulations) repository which is an integrated framework for computational science using density functional theory, classical force-field/molecular dynamics and machine-learning. The NIST-JARVIS official website is: https://jarvis.nist.gov . This project is a part of the Materials Genome Initiative (MGI) at NIST (https://mgi.nist.gov/). 

For more details, checkout our latest article:  `The joint automated repository for various integrated simulations (JARVIS) for data-driven materials design <https://www.nature.com/articles/s41524-020-00440-1>`__ and `YouTube videos <https://www.youtube.com/watch?v=P0ZcHXOC6W0&feature=emb_title&ab_channel=JARVIS-repository>`__ 



Capabilities
---------------------

- **Software workflow tasks for preprcessing, executing and post-processing**:  VASP, Quantum Espresso, Wien2k BoltzTrap, Wannier90, LAMMPS, Scikit-learn, TensorFlow, LightGBM, Qiskit, Tequila, Pennylane, DGL, PyTorch.

- **Several examples**: Notebooks and test scripts to explain the package.

- **Several analysis tools**: Atomic structure, Electronic structure, Spacegroup, Diffraction, 2D materials and other vdW bonded systems, Mechanical, Optoelectronic, Topological, Solar-cell, Thermoelectric, Piezoelectric, Dielectric, STM, Phonon, Dark matter, Wannier tight binding models, Point defects, Heterostructures, Magnetic ordering, Images, Spectrum etc.

- **Database upload and download**: Download JARVIS databases such as JARVIS-DFT, FF, ML, WannierTB, Solar, STM and also external databases such as Materials project, OQMD, AFLOW etc.

- **Access raw input/output files**: Download input/ouput files for JARVIS-databases to enhance reproducibility.

- **Train machine learning models**: Use different descriptors, graphs and datasets for training machine learning models.

- **HPC clusters**: Torque/PBS and SLURM.


Installation
---------------------
Using pip
^^^^^^^^^
>>> pip install -U jarvis-tools

or

Using conda
^^^^^^^^^^^^
First create a conda environment: Install miniconda environment from https://conda.io/miniconda.html Based on your system requirements, you'll get a file something like 'Miniconda3-latest-XYZ'.

>>> bash Miniconda3-latest-Linux-x86_64.sh (for linux)
>>> bash Miniconda3-latest-MacOSX-x86_64.sh (for Mac)

Now let's create a conda environment and install jarvis-tools from conda-forge:

>>> conda create --name my_jarvis python=3.8
>>> source activate my_jarvis
>>> conda install -c conda-forge jarvis-tools

Please make sure to use python>3.7.


Example function
---------------------
Create Atoms object:

>>> from jarvis.core.atoms import Atoms
>>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
>>> coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
>>> elements = ["Si", "Si"]
>>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
>>> density = round(Si.density,2)
>>> print (density)
2.33

Obtain JARVIS-DFT 3D dataset with various materials and their properties

>>> from jarvis.db.figshare import data
>>> dft_3d = data(dataset='dft_3d')
>>> print (len(dft_3d))
48527

Write to POSCAR files to visualize/analyze in VESTA or other packages

>>> from jarvis.io.vasp.inputs import Poscar
>>> for i in dft_3d:
...     atoms = Atoms.from_dict(i['atoms'])
...     poscar = Poscar(atoms)
...     jid = i['jid']
...     filename = 'POSCAR-'+jid+'.vasp'
...     poscar.write_file(filename)

Get JARVIS-DFT 2D dataset

>>> dft_2d = data(dataset='dft_2d')
>>> print (len(dft_2d))
1070
>>> for i in dft_2d:
...     atoms = Atoms.from_dict(i['atoms'])
...     poscar = Poscar(atoms)
...     jid = i['jid']
...     filename = 'POSCAR-'+jid+'.vasp'
...     poscar.write_file(filename)



Example to parse DOS data from JARVIS-DFT XML webpages


>>> from jarvis.db.webpages import Webpage
>>> from jarvis.core.spectrum import Spectrum
>>> import numpy as np
>>> new_dist=np.arange(-5, 10, 0.05)
>>> all_atoms = []
>>> all_dos_up = []
>>> all_jids = []
>>> for ii,i in enumerate(dft_3d):
      all_jids.append(i['jid'])
...   try:
...     w = Webpage(jid=i['jid'])
...     edos_data = w.get_dft_electron_dos()
...     ens = np.array(edos_data['edos_energies'].strip("'").split(','),dtype='float')
...     tot_dos_up = np.array(edos_data['total_edos_up'].strip("'").split(','),dtype='float')
...     s = Spectrum(x=ens,y=tot_dos_up)
...     interp = s.get_interpolated_values(new_dist=new_dist)
...     atoms=Atoms.from_dict(i['atoms'])
...     ase_atoms=atoms.ase_converter()
...     all_dos_up.append(interp)
...     all_atoms.append(atoms)
...     all_jids.append(i['jid'])
...     filename=i['jid']+'.cif'
...     atoms.write_cif(filename)
...     break
...   except Exception as exp :
...     print (exp,i['jid'])
...     pass



Find more examples at

      
      1) https://jarvis-tools.readthedocs.io/
      
      2) https://github.com/JARVIS-Materials-Design/jarvis-tools-notebooks
      
      3) https://github.com/usnistgov/jarvis/tree/master/jarvis/tests/testfiles
      
Citing
---------------------

Please cite the following if you happen to use JARVIS-Tools for a publication.

https://www.nature.com/articles/s41524-020-00440-1

  @article{choudhary2020joint,
    title={The joint automated repository for various integrated simulations (JARVIS) for data-driven materials design},
    author={Choudhary, Kamal and Garrity, Kevin F and Reid, Andrew CE and DeCost, Brian and Biacchi, Adam J and Walker, Angela R Hight and Trautt, Zachary and Hattrick-Simpers, Jason and Kusne, A Gilad and Centrone, Andrea and others},
    journal={npj Computational Materials},
    volume={6},
    number={1},
    pages={1--13},
    year={2020},
    publisher={Nature Publishing Group}
  }

      

Module details
--------------------

* :ref:`modindex`
* :ref:`genindex`

Correspondence
--------------------

Please report bugs as Github issues (https://github.com/usnistgov/jarvis/issues) or email to kamal.choudhary@nist.gov.

Funding support
--------------------

NIST-MGI (https://www.nist.gov/mgi).

Code of conduct
--------------------

Please see `Code of conduct <https://github.com/usnistgov/jarvis/blob/master/CODE_OF_CONDUCT.md>`__

