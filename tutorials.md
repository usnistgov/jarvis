# Tutorials

## Overview

This page introduces several tools and dataset in JARVIS though examples.
In addition to the following examples, the [Colab notebooks](https://github.com/JARVIS-Materials-Design/jarvis-tools-notebooks)
and [pytests modules](https://github.com/usnistgov/jarvis/tree/master/jarvis/tests/testfiles) can be helpful to 
get use JARVIS-DB and JARVIS-Tools.


## DFT

### Running calculations
The JARVIS-Tools curretly allows to run DFT calculations with VASP and QE software. The JARVIS-DFT is mainly based on VASP software but 
sooon there would be datasets with QE asl well. We can create for example a VaspJob with the help of
atomic structure, input parameters, pseudopotential, k-points information. Similar to many other modules
VaspJob allows to 'ToDict' and 'FromDict' methods to store or load a complete job, which is very useful in
scaling up VASP related calculations and enhancing reproducibilty.
Make sure JARVIS_VASP_PSP_DIR is declared as a PATH to VASP pseudopotential directory.
The input file generation and output file parsing modules for VASP can be found 
in jarvis.io.vasp.inputs and jarvis.io.vasp.outputs modules.
We start by setting up and submitting a single VaspJob:
``` python hl_lines="3"
from jarvis.tasks.vasp.vasp import VaspJob, write_vaspjob
from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
import os

# Load/build crystal structure
# mat = Poscar.from_file('POSCAR')
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
elements = ["Si", "Si"]
box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
mat = Poscar(atoms)
mat.comment = "Silicon"

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
# export JARVIS_VASP_PSP_DIR = 'PATH_TO_YOUR_PSP'
pot = Potcar(elements=mat.atoms.elements)

# Build Kpoints info
kp = Kpoints3D().automatic_length_mesh(
    lattice_mat=mat.atoms.lattice_mat, length=20
)

vasp_cmd = "/users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_std"
copy_files = ["/users/knc6/bin/vdw_kernel.bindat"]
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


``` 

The job.py can now be run on a cluster or on a PC as a python script.
For running this job on a PBS cluster,

``` python hl_lines="3"
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

 
``` 
Currently, JARVIS-Tools can be used to submit job with SLURM and PBS clusters only.
For high-throughput automated submissions one can use pre-build JobFactory module
that allows automatic calculations for a series of properties.
``` python hl_lines="3"

from jarvis.tasks.vasp.vasp import VaspJob,write_jobfact_optb88vdw
from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
import os
from jarvis.tasks.vasp.vasp import JobFactory,write_jobfact_optb88vdw

# Load/build crystal structure
# mat = Poscar.from_file('POSCAR')
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
elements = ["Si", "Si"]
box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
mat = Poscar(atoms)
mat.comment = 'bulk@Silicon'

vasp_cmd = "/users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_std"
copy_files = ['/users/knc6/bin/vdw_kernel.bindat']
jobname = "MAIN-RELAX@JVASP-1002"
job = JobFactory(
    vasp_cmd=vasp_cmd,
    poscar=mat,
    copy_files=copy_files,
)

dumpjson(data=job.to_dict(), filename='job_fact.json')
write_jobfact_optb88vdw(pyname="job_fact.py", job_json="job_fact.json")

``` 
We can now submit the JobFactory as:
``` python hl_lines="3"
from jarvis.tasks.queue_jobs import Queue
import os
submit_cmd=["qsub", "submit_job"]
# Example job commands, need to change based on your cluster
job_line = "source activate my_jarvis \n"+"python job_fact.py"
name = "TestJob"
directory = os.getcwd()
Queue.pbs(
    job_line=job_line,
    jobname=name,
    directory=directory,
    submit_cmd=submit_cmd,
)


``` 

This script first converges K-pints and plane-wave cutoff, then using the converged paramerters
optimizes the simulation cell, and then runs several jobs for calculating properties such as 
bandstructure on high-symmetry k-points, elastic constanats, optoelectronic properties etc.

Next, let's see how to run High-throughput jobs for multiple JARVIS-IDs:

``` python hl_lines="3"
# Complete workflow for running high-throughput calculations
from jarvis.tasks.vasp.vasp import VaspJob, write_jobfact_optb88vdw, JobFactory
from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
import os
from jarvis.db.figshare import get_jid_data

#aoth to executable
vasp_cmd = "mpirun /home/knc6/Software/vasp.5.4.1/bin/vasp_std"
# Needed for vdW functionals only, else keep it an empty array
copy_files = ["/home/knc6/bin/vdw_kernel.bindat"]
# submit_cmd = ["qsub", "submit_job"] # For Torque clusters
submit_cmd = ["sbatch", "submit_job"]
# List of JARVIS-IDs for which structure would be fetched and subjected to high-throughput
jids = ["JVASP-1002", "JVASP-1067"]

for jid in jids:
    d = get_jid_data(jid=jid, dataset="dft_3d")
    # Make Atoms class from python dictinary object
    atoms = Atoms.from_dict(d["atoms"])
    mat = Poscar(atoms)
    mat.comment = "bulk@" + str(jid)
    cwd_home = os.getcwd()
    dir_name = d["jid"] + "_" + str("PBEBO")
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)
    job = JobFactory(
        vasp_cmd=vasp_cmd,
        poscar=mat,
        copy_files=copy_files,
    )
    dumpjson(data=job.to_dict(), filename="job_fact.json")
    write_jobfact_optb88vdw(pyname="job_fact.py", job_json="job_fact.json")

    # Example job commands, needs to be changed based on your cluster
    job_line = (
        "source activate my_jarvis \n"
        + "module load intel/2015 openmpi/1.10.2/intel-15 \n"
        + "python job_fact.py"
    )
    name = jid
    directory = os.getcwd()
    Queue.slurm(
        job_line=job_line,
        jobname=name,
        directory=directory,
        submit_cmd=submit_cmd,
    )
    os.chdir(cwd_home)
    """
    # For Torque clusters
    Queue.pbs(
        job_line=job_line,
        jobname=name,
        directory=directory,
        submit_cmd=submit_cmd,
    )
    os.chdir(cwd_home)
    """
```
In the above examples, use Queue.slurm if you want to use SLURM instead of TORQUE/PBS submission.
A complete example of such run is available at: [VASP example](https://github.com/usnistgov/jarvis/tree/master/jarvis/examples/vasp)

### Post-processing and plotting
There are a variety of post-processing analysis and plotting that can be done on the output data.
A common example would be plotting electronic density of states and bandstructure as follows:
``` python hl_lines="3"
from jarvis.io.vasp.outputs import Vasprun
# %matplotlib inline
import matplotlib.pyplot as plt
plt.switch_backend('agg')
vrun = Vasprun('vasprun.xml')
energies, spin_up, spin_dn=vrun.total_dos
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.plot(energies,spin_up,label='Spin-up')
plt.plot(energies,spin_dn,label='Spin-down')
plt.xlabel('Energy(E-Ef)')
plt.ylabel('DOS(arb.unit)')
plt.xlim(-4,4)
plt.legend()

vrun.get_bandstructure(kpoints_file_path='KPOINTS')
``` 
There are many other modules available such as: scanning tunneling microscopy images, solar-cell efficiency, topological spin-
orbit spillage, transport properties, phonons, infrared intensities, and its continously expanding.


### Developing database

After generating the results, we can store the metadata in JARVIS-API. Please request an account if you haven't made it yet.
Following the calculation protocal mentioned above, the generated files can be converted to an XML datafile which with the help of
an XSD schema can be converted to nice-looking HTML files with the help of XSLT programming and a bit of javascript.

We alreay provide modules to convert the calculation informato to XML and module to upload data. An example is give below:


``` python hl_lines="3"
from jarvis.db.vasp_to_xml import VaspToApiXmlSchema 
from jarvis.db.restapi import Api
folder="/home/users/knc6/Software/jarvis/jarvis/examples/vasp/SiOptB88vdW"
filename = "JVASP-1002.xml"
VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)
a = Api(user_info_file_path="/users/knc6/myinfo")
# First line should be your username and secondline your password
tid="5f626925ece4b00035e5277f"
# Find latest template ID with title "jarvisdft"  at the bottom of the page 
# https://jarvis.nist.gov/rest/template-version-manager/global
a.upload_xml_file(filename='JVASP-1067.xml',template_id=tid)
``` 

## FF

Molecular dynamics/classical force-field calculations can be carried out with LAMMPS software.
An example for running LAMMPS is given below. Here, a LammpsJob module is defined with the help of 
atoms, pair-style, coefficient, and template file (*.mod file) to control the calculations.

### Running calculations
``` python hl_lines="3"
from jarvis.tasks.lammps.lammps import LammpsJob, JobFactory
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from jarvis.analysis.structure.spacegroup import Spacegroup3D


# atoms = Atoms.from_poscar('POSCAR')
# Get Aluminum FCC
tmp_dict = get_jid_data(jid="JVASP-816", dataset="dft_3d")["atoms"]
atoms = Atoms.from_dict(tmp_dict)
# Get conventional cell
spg = Spacegroup3D(atoms)
cvn_atoms = spg.conventional_standard_structure
ff = "/users/knc6/Software/LAMMPS/lammps-master/potentials/Al_zhou.eam.alloy"
mod = "/users/knc6/Software/Devs/jarvis/jarvis/tasks/lammps/templates/inelast.mod"
cmd = "/users/knc6/Software/LAMMPS/lammps-master/src/lmp_serial<in.main>out"
parameters = {
    "pair_style": "eam/alloy",
    "pair_coeff": ff,
    "atom_style": "charge",
    "control_file": mod,
}
# Test LammpsJob
lmp = LammpsJob(
    atoms=cvn_atoms, parameters=parameters, lammps_cmd=cmd, jobname="Test"
).runjob()

# Test high-throughput
job_fact = JobFactory(pair_style="eam/alloy", name="my_first_lammps_run")
job_fact.all_props_eam_alloy(atoms=cvn_atoms, ff_path=ff, lammps_cmd=cmd)
```
The above lines can be written on a file such as job.py and can be run on a PC or a cluster with the
from jarvis.tasks.queue_jobs import Queue module.


### Post-processing and plotting

Important quantities such asn total energy, forces etc. can be obtained with the help of
jarvis.io.lammps.outputs module. 

``` python hl_lines="3"
from jarvis.io.lammps.outputs import parse_material_calculation_folder
folder = '/home/users/knc6/Software/jarvis/jarvis/examples/lammps/Aleam'
data = parse_material_calculation_folder(folder)
``` 

### Developing database

The calculation data can now be converted into XML files as follows:

``` python hl_lines="3"
from jarvis.db.lammps_to_xml import write_xml
write_xml(data=data,filename='lmp.xml')
``` 

## ML/AI

Currently JARVIS-ML allows prediction of material properties with machine learning. The materials information
is converted into descriptors using Classical Force-field Inspired Descriptors (CFID) or Coulomb materix.
Other descriptors and graph based predictions would be available soon also.
For a series of atomistic structures, we can convert them into CFID, which act as input matrix.

### Running calculations
Suppose we have 40000 materials, and we get 1557 descriptor for each material, we wull have a 
40000x1557 matrix. Let's call this matrix as 'x' or input matrix.
Next, we can get target ('y') data either from DFT, FF calculations or experiments.
For example, we can choose formation energies of 40000 materials in the JARVIS-DFT as the dtarget data
giving 40000x1 matrix.

Now, we can use a ML/AI algorithm to establish statistical relation between the x and y data.
Once trained we get a trained model, which can be stored in say pickle or joblib format.

For a new material now, it can be converted into CFID i.e. 1x1557 matrix which when fed to the model
will give 1x1 prediction hence the ML prediction. We can use a range of ML algorithms such as 
linear regression, decision trees, Gaussian processes etc. We find with CFID descriptors, gradient boosting
decision trees (especially in LightGBM) gives one of the most accurate results. 
We provide tools to run with major ML packages such as scikit-learn, tensorflow, pytorch, lightgbm etc.

``` python hl_lines="3"
# An example of JARVIS-ML training
from jarvis.ai.pkgs.utils import get_ml_data
from jarvis.ai.pkgs.utils import regr_scores
X,y,jid=get_ml_data()
#Formation energy for 3D materials, you can choose other properties/dataset as well
import lightgbm as lgb
from sklearn.model_selection import train_test_split
lgbm = lgb.LGBMRegressor(device= 'gpu',n_estimators= 1170,learning_rate= 0.15375236057119931,num_leaves= 273)       
X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(X, y, jid, random_state=1, test_size=.1)
lgbm.fit(X_train,y_train)
pred = lgbm.predict(X_test)
reg_sc = regr_scores(y_test, pred)
print (reg_sc['mae'])

```

Another full example for regression
``` python hl_lines="3"
from jarvis.ai.pkgs.lgbm.regression import parameters_dict
from scipy.stats import median_absolute_deviation as mad
from jarvis.ai.pkgs.utils import get_ml_data
import lightgbm as lgb
import numpy as np
from sklearn.model_selection import train_test_split
from jarvis.ai.pkgs.utils import regr_scores
import joblib

params = parameters_dict()

print(params)


mem = []
for i, j in params.items():
    name = str(i) + ".pkl"
    print(i)
    print(name)
    X, y, jid = get_ml_data(dataset="cfid_3d", ml_property=i)
    lgbm = lgb.LGBMRegressor(
        n_estimators=j["n_estimators"],
        learning_rate=j["learning_rate"],
        num_leaves=j["num_leaves"],
    )
    if "eps" in i:  # fit refractive index, not dielectric constant
        y = np.sqrt(y)
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        X, y, jid, random_state=1, test_size=0.1
    )
    lgbm.fit(X_train, y_train)
    pred = lgbm.predict(X_test)
    joblib.dump(lgbm, name)
    reg_sc = regr_scores(y_test, pred)
    mae = reg_sc["mae"]  # mean absolute error
    madev = mad(y)  # mean absolute deviation
    # mae_over_madev=float(mae)/float(madev)
    mem.append([i, len(X), mae, madev])
    print("Property,Length, MAE,MAD", i, len(X), mae, madev)
    print()
    print()
```
### Post-processing and plotting

We can analyze basic ML metrics such as mean-absolute erros (MAE),  RMSE, R2 etc. for regression and
ROC AUC, F1 score etc. for classification models. 

``` python hl_lines="3"
from jarvis.ai.pkgs.sklearn.regression import regression
from jarvis.ai.pkgs.sklearn.classification import classification
from jarvis.ai.pkgs.sklearn.hyper_params import (
    classification_regression_params,
)
from jarvis.ai.pkgs.utils import get_ml_data, binary_class_dat
from jarvis.ai.pkgs.lgbm.regression import regression as l_regression
from jarvis.ai.pkgs.lgbm.regression import parameters_dict as l_params
from jarvis.ai.pkgs.lgbm.classification import (
    classification as l_classification,
)
from jarvis.ai.descriptors.cfid import feat_names
from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt
property = "exfoliation_energy"
X, Y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)

# Regression
params = l_params()[property]
names = feat_names()
info = l_regression(
    X=X,
    Y=Y,
    jid=jid,
    config=params,
    feat_names=names,
    plot=True,
    save_model=True,
)

# Classification
property = "exfoliation_energy"
tol = 100
models = [LGBMClassifier(n_estimators=10, num_leaves=2)]
info = l_classification(
    X=X,
    Y=Y,
    plot=True,
    models=models,
    preprocess=True,
    save_model=True,
    tol=tol,
)
```
### Developing database
The ML models in terms of joblib, pickle parameters, the descriptors sets and ML predicted properties
can be stored in teh JARVIS-API.
The descriptrs and the predicted properties can be stored as XML files , while the joblib and pickle files
can be saved a blob binary files.



