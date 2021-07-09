Tutorials
=============

.. _customise-templates:


How to analyze an atomic structure
------------------------------------------------------------

How to get RDF, ADF, DDF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get XRD paterns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to setup/analyze DFT calculations using VASP
------------------------------------------------
The Vienna Ab initio Simulation Package, better known as VASP, is a package for performing ab initio quantum mechanical calculations using either Vanderbilt pseudopotentials, or the projector augmented wave method, and a plane wave basis set. Manual for VASP is available at: https://www.vasp.at/wiki/index.php/The_VASP_Manual .

Running a VASP calculation requires the following files: ``INCAR``, ``POSCAR``, ``KPOINTS``, ``POTCAR`` as well as additional files such as ``vdw_kernel.bindat`` for specific types of calculations. While setting up calculations for one or a few systems/setups should be straight forward, setting up calculations for thousands of materials and most importantly making a database out of all those calculations require automated calculations script collections such as JARVIS-Tools. 

Gievn an atomic structure in 1) ``jarvis.core.Atoms`` format, JARVIS-Tools 2) prepares input files such as ``INCAR`` etc. as mentioned above and 3) submits the calculations to your queuing system such as SLURM/PBS using ``jarvis.tasks.vasp`` and ``jarvis.tasks.queue_jobs``. After a calculations get completed, 4) automated analysis can be carried out and plots and webpages are generated. The input file generation and output file parsing modules for VASP can be found in ``jarvis.io.vasp.inputs`` and ``jarvis.io.vasp.outputs`` modules. The automated analyis and XML generation for webpages can be found in ``jarvis.db.vasp_to_xml`` module. After the xml page creation they are converted using html using XSLT scripts. 

Additionally, a JSON file is created with metadata from all the XML pages for thousands of materials to easily use in data-analytics/machine learning applications.The JARVIS-DFT (https://jarvis.nist.gov/jarvisdft/) database primarily uses such a workflow.
Make sure `VASP_PSP_DIR` is declared as a PATH to VASP pseudopotential directory i.e. 

.. highlight:: bash

::

   $ export VASP_PSP_DIR=YOUR_PATH_TO_PSUEDOPTENTIALS

in your ~/.bashrc file.

How to setup a single calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We start by setting up and submitting a single VaspJob:

.. code-block:: python


   from jarvis.tasks.vasp.vasp import VaspJob, write_vaspjob
   from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
   from jarvis.db.jsonutils import dumpjson
   from jarvis.core.atoms import Atoms
   from jarvis.core.kpoints import Kpoints3D
   from jarvis.tasks.queue_jobs import Queue
   import os

   # Load/build crystal structure
   mat = Poscar.from_file('POSCAR')
   # coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
   # elements = ["Si", "Si"]
   # box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
   # atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
   # mat = Poscar(atoms)
   # mat.comment = "Silicon"

   # Build INCAR file
   data = dict(
       PREC="Accurate",
       ISMEAR=0,
       SIGMA=0.01,
       IBRION=2,
       ISIF=3,
       GGA="BO",
       PARAM1=0.1833333333,
       PARAM2=0.2200000000,
       LUSE_VDW=".TRUE.",
       AGGAC=0.0000,
       EDIFF="1E-7",
       EDIFFG="-1E-3",
       NELM=400,
       ISPIN=2,
       LCHARG=".FALSE.",
       LVTOT=".FALSE.",
       LVHAR=".FALSE.",
       LWAVE=".FALSE.",
   )
   inc = Incar(data)
   # Build POTCAR info
   # export VASP_PSP_DIR = 'PATH_TO_YOUR_PSP'
   pot = Potcar.from_atoms(mat.atoms)
   #pot = Potcar(elements=mat.atoms.elements)

   # Build Kpoints info
   kp = Kpoints3D().automatic_length_mesh(
       lattice_mat=mat.atoms.lattice_mat, length=20
   )

   vasp_cmd = "PATH_TO vasp_std"
   copy_files = ["PATH_TO vdw_kernel.bindat"]
   jobname = "MAIN-RELAX@JVASP-1002"
   job = VaspJob(
       poscar=mat,
       incar=inc,
       potcar=pot,
       kpoints=kp,
       vasp_cmd=vasp_cmd,
       copy_files=copy_files,
       jobname=jobname,
   )

   dumpjson(data=job.to_dict(), filename="job.json")
   write_vaspjob(pyname="job.py", job_json="job.json")

The job.py can now be run on a cluster or on a PC as a python script. For running this job on a PBS cluster,

.. code-block:: python


   submit_cmd = ["qsub", "submit_job"]
   # Example job commands, need to change based on your cluster
   job_line = (
       "source activate my_jarvis \n"
       + "python job.py"
   )
   name = "TestJob"
   directory = os.getcwd()
   Queue.pbs(
       job_line=job_line,
       jobname=name,
       directory=directory,
       submit_cmd=submit_cmd,
       )





How to setup high-throughput calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, JARVIS-Tools can be used to submit job with SLURM and PBS clusters only. For high-throughput automated submissions one can use pre-build JobFactory module that allows automatic calculations for a series of properties.



.. code-block:: python


   # List of materials to run high-throughput calculations on
   ids = ['POSCAR-1.vasp','POSCAR-2.vasp','POSCAR-3.vasp']

   from jarvis.tasks.vasp.vasp import (
       JobFactory,
       VaspJob,
       GenericIncars,
       write_jobfact,
   )
   from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
   from jarvis.db.jsonutils import dumpjson
   from jarvis.db.figshare import data
   from jarvis.core.atoms import Atoms
   from jarvis.tasks.queue_jobs import Queue
   import os
   vasp_cmd = "mpirun PATH_TO vasp_std"
   copy_files = ["PATH_TO vdw_kernel.bindat"]
   submit_cmd = ["qsub", "submit_job"]

   # For slurm
   # submit_cmd = ["sbatch", "submit_job"]

   steps = [
       "ENCUT",
       "KPLEN",
       "RELAX",
       "BANDSTRUCT",
       "OPTICS",
       "MBJOPTICS",
       "ELASTIC",
   ]
   incs = GenericIncars().optb88vdw().incar.to_dict()

   for id in ids:
       mat = Poscar.from_file(id)
       cwd_home = os.getcwd()
       dir_name = id.split('.vasp')[0] + "_" + str("PBEBO")
       if not os.path.exists(dir_name):
           os.makedirs(dir_name)
       os.chdir(dir_name)
       job = JobFactory(
           vasp_cmd=vasp_cmd,
           poscar=mat,
           steps=steps,
           copy_files=copy_files,
           use_incar_dict=incs,
       )

       dumpjson(data=job.to_dict(), filename="job_fact.json")
       write_jobfact(
           pyname="job_fact.py",
           job_json="job_fact.json",
           input_arg="v.step_flow()",
       )

       # Example job commands, need to change based on your cluster
       job_line = (
           "source activate my_jarvis \n"
           + "python job_fact.py"
       )
       name = id
       directory = os.getcwd()
       Queue.pbs(
           job_line=job_line,
           jobname=name,
           #partition="",
           walltime="24:00:00",
           #account="",
           cores=12,
           directory=directory,
           submit_cmd=submit_cmd,
       )
       os.chdir(cwd_home)
       """
       # For Slurm clusters
       Queue.slurm(
           job_line=job_line,
           jobname=name,
           directory=directory,
           submit_cmd=submit_cmd,
       )
       os.chdir(cwd_home)
       """


We provide modules to convert the calculation informato to ``XML`` which can be converted to ``HTML`` using ``XSLT``. An example is give below:

.. code-block:: python

   from jarvis.db.vasp_to_xml import VaspToApiXmlSchema
   from jarvis.db.restapi import Api
   folder="jarvis/jarvis/examples/vasp/SiOptB88vdW"
   filename = "JVASP-1002.xml"
   VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)


How to plot electronic bandstructure and DOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you use the workflow used above, the density of states plot can be obtained using thr ``vasprun.xml`` file in MAIN-RELAX folder while the band-structure plot is obtained using ``vasprun.xml`` in MAIN-BAND folder.

.. code-block:: python


   from jarvis.io.vasp.outputs import Vasprun
   vrun = Vasprun('vasprun.xml')
   %matplotlib inline
   import matplotlib.pyplot as plt
   plt.rcParams.update({'font.size': 22})

   # Bandstructure plot
   vrun.get_bandstructure(kpoints_file_path='KPOINTS')

   # DOS plot
   energies, spin_up, spin_dn=vrun.total_dos
   plt.rcParams.update({'font.size': 22})
   plt.plot(energies,spin_up,label='Spin-up')
   plt.plot(energies,spin_dn,label='Spin-down')
   plt.xlabel('Energy(E-Ef)')
   plt.ylabel('DOS(arb.unit)')
   plt.xlim(-4,4)
   plt.legend()


How to obtain elastic constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to plot generate an STM/STEM image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to plot generate a dielectric function spectra and solar eff.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to generate/use electronic Wannier tight binding model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to generate Fermi-surfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to run BoltzTrap for transport properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to make heterostructures/interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get IR/Raman spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get piezoelectic/dielecrric/BEC constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get electric field gradients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get work-function of a surface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to get exfoliation energy of a 2D material
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to run/analyze MD static/dynamic calculation using LAMMPS
-------------------------------------------------------------

How to setup a single calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to setup high-throughput calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to setup computer-cluser details
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


How to run/analyze DFT static calculation using Quantum espresso
-----------------------------------------------------------------

How to setup a single calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to setup high-throughput calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to setup computer-cluser details
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



How to traing JARVIS-CFID ML models using sklearn/lightgbm
----------------------------------------------------------

How to train regression model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to train classification model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


How to traing JARVIS-ALIGNN ML models using PyTorch
-----------------------------------------------------

How to train regression model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to train classification model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


How to use quantum computation algorithms using Qiskit/Tequila/Pennylane
------------------------------------------------------------------------

How to generate circuit model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to run cals. on simulators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How to run cals. on actual quantum computers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
