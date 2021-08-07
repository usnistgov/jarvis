Databases
=============

.. _customise-templates:

====================     =========================   ======================================================= 
Database name            Number of data-points       Description                                             
====================     =========================   ======================================================= 
``dft_3d``               48527                       Various 3D materials properties in JARVIS-DFT database  
                                                     computed with OptB88vdW and TBmBJ methods             
``dft_2d``               1079                        Various 2D materials properties in JARVIS-DFT database  
                                                     computed with OptB88vdW                                
``stm``                  1132                        2D materials STM images in JARVIS-STM database  
``wtbh_electron``        1440                        3D and 2D materials Wannier tight-binding Hamiltonian
                                                     dtaabase for electrons with spin-orbit coupling
                                                     in JARVIS-WTB (Keyword: 'WANN')
``wtbh_phonon``          15502                       3D and 2D materials Wannier tight-binding Hamiltonian
                                                     for phonons at Gamma with finite difference 
                                                     (Keyword:FD-ELAST)
``jff``                  2538                        Various 3D materials properties in JARVIS-FF database   
                                                     computed with several force-fields                     
``edos_pdos``            48469                       Normalized electron and phonon density of states with 
                                                     interpolated values and fixed number of bins
``megnet``               69239                       Formation energy and bandgaps of 3D materials properties
                                                     in Materials project database as on 2018, used in megnet
``twod_matpd``           6351                        Formation energy and bandgaps of 2D materials properties
                                                     in 2DMatPedia database
``c2db``                 3514                        Various properties in C2DB database
``polymer_genome``       1073                        Electronic bandgap and diecltric constants of crystall
                                                     ine polymer in polymer genome database
``qm9_std_jctc``         130829                      Various properties of molecules in QM9 database
``cod``                  431778                      Atomic structures from crystallographic open database
``oqmd_3d_no_cfid``      817636                      Formation energies and bandgaps of 3D materials 
                                                     from OQMD database
``omdb``                 12500                       Bandgaps  for organic polymers in OMDB database
``hopv``                 4855                        Various properties of molecules in HOPV15 dataset 
``pdbbind``              11189                       Bio-molecular complexes database from PDBBind v2015
``qmof``                 18321                       Bandgaps and total energies of metal organic frameowrks
                                                     in QMOF database
``raw_files``            144895                      Figshare links to download raw calculations VASP files
                                                     from JARVIS-DFT
====================     =========================   ======================================================= 

All these datasets can be obtained using jarvis-tools as follows, exception to ``stm``, ``wtbh_electron``, ``wtbh_phonon``
which have their own modules in ``jarvis.db.figshare``:

.. code-block:: python

                from jarvis.db.figshare import data
                d = data('dft_3d') #choose a name of dataset from above
                # See available keys
                print (d[0].keys())
                # Dataset size
                print (len(d)

                # Visualize an atoms object
                from jarvis.core.atoms import Atoms
                a = Atoms.from_dict(d[0]['atoms'])
                #You can visualize this in VESTA or other similar packages
                print (a)

                # If pandas framework needed
                import pandas as pd
                df = pd.DataFrame(d)
                print (df)

JARVIS-DFT
------------------------------------------------


JARVIS-Formation energy and bandgap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-2D Exfoliation energies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-MetaGGA (dielectric function and SLME, solar cells)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-STM and STEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-WannierTB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-Elastic constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-Topological materials (Spin-orbit Spillage)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-DFPT (Piezoelectric, IR, Raman, dielectric, BEC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-BoltzTrap (Thermoelectrics coeff, eff. mass)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-Magnetic moments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-DFPT (Piezoelectric, IR, dielectric)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-EFG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-PBE0 and HSE06
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-Heterostructure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-EDOS-PDOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-Kpoint and cut-off
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-FF
-------------------------------------------------------------

Energetics
^^^^^^^^^^^^

Elastic constants
^^^^^^^^^^^^

Vacancy formation energy
^^^^^^^^^^^^

Surface energy and Wulff-plots
^^^^^^^^^^^^

Phonon DOS
^^^^^^^^^^^^

JARVIS-RAW Files
-------------------------------------------------------------

JARVIS-DFT structure relaxation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-DFT Elastic constants/finite difference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-WannierTB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JARVIS-STM and STEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

External datasets used for ML training
-------------------------------------------------------------

Materials project dataset 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

QM9 dataset 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OQMD dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AFLOW dataset 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Polymer genome dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

COD dataset 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OMDB dataset 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

QMOF dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C2DB dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

HPOV dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
