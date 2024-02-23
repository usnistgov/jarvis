# Tutorials

## How to analyze an atomic structure

Atomic structure act as an input to multiple simulations such as for
density functional theory, molecular dyanmics, Monte Carlo, atomistic
graph neural network etc. So, we provide a very bried introduction to
the atomic structure here. For more general information, refer to
solid-state physics or introduction to materials-science books.

An atomic structure can consist of atomic element types, corresponding
xyz coordinates in space (either in real or reciprocal space) and
lattice matrix used in setting periodic boundary conditions.

An example of constructing an atomic structure class using
`jarvis.core.Atoms` is given below. After creating the Atoms class, we
can simply print it and visualize the POSCAR format file in a software
such as VESTA. While the examples below use Silicon elemental crystal
creation and analysis, it can be used for multi-component systems as
well.

``` python
from jarvis.core.atoms import Atoms
box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
elements = ["Si", "Si"]
Si = Atoms(lattice_mat=box, coords=coords, elements=elements, cartesian=False)
print (Si) # To visualize 
Si.write_poscar('POSCAR.vasp')
Si.write_cif('POSCAR.vasp')
```

The <span class="title-ref">Atoms</span> class here is created from the
raw data, but it can also be read from different file formats such as:
<span class="title-ref">'.cif', 'POSCAR', '.xyz', '.pdb', '.sdf',
'.mol2'</span> etc. The Atoms class can also be written to files in
formats such as POSCAR/.cif etc.

Note that for molecular systems, we use a large vaccum padding (say 50
Angstrom in each direction) and set lattice_mat accordingly, e.g.
lattice_mat = \[\[50,0,0\],\[0,50,0\],\[0,0,50\]\]. Similarly, for free
surfaces we set high vaccum in one of the crystallographic directions
(say z) by giving a large z-comonent in the lattice matrix while keeping
the x, y comonents intact.

``` python
my_atoms = Atoms.from_poscar('POSCAR')
my_atoms.write_poscar('MyPOSCAR')
```

Once this Atoms class is created, several imprtant information can be
obtained such as:

``` python
print ('volume',Si.volume)
print ('density in g/cm3', Si.density)
print ('composition as dictionary', Si.composition)
print ('Chemical formula', Si.composition.reduced_formula)
print ('Spacegroup info', Si.spacegroup())
print ('lattice-parameters', Si.lattice.abc, Si.lattice.angles)
print ('packing fraction',Si.packing_fraction)
print ('number of atoms',Si.num_atoms)
print ('Center of mass', Si.get_center_of_mass())
print ('Atomic number list', Si.atomic_numbers)
```

For creating/accessing dataset(s), we use `Atoms.from_dict()` and
`Atoms.to_dict()` methods:

``` python
d = Si.to_dict()
new_atoms = Atoms.from_dict(d)
```

The <span class="title-ref">jarvis.core.Atoms</span> object can be
converted back and forth to other simulation toolsets such as Pymatgen
and ASE if insyalled, as follows

``` python
pmg_struct = Si.pymatgen_converter()
ase_atoms = Si.ase_converter()
```

In order to make supercell, the following example can be used:

``` python
supercell_1 = Si.make_supercell([2,2,2])
supercell_2 = Si.make_supercell_matrix([[2,0,0],[0,2,0],[0,0,2]])
supercell_1.density == supercell_2.density
```

### How to get RDF, ADF, DDF

Nearest-neighbor analysis one of the most important tools in atomistic
simulations. Quantities such as radial (RDF), angle (ADF) and dihedral
(DDF) distribution functions can be obtained using
<span class="title-ref">jarvis.analysis.structure.neighbors.NeighborsAnalysis</span>
class as shown in the following example using the Si Atoms class
obtained above. Different cut-off parameters for angle and sihedral
distribution are used to narrow down the number of neighbors. For
details, please look into respective modules.

``` python
nb = NeighborsAnalysis(Si)
bins_rdf, rdf, nbs = nb.get_rdf() #Global Radial distribution function
adfa, bins_a = nb.ang_dist_first() #Angular distribution function upto first neighbor
adfb, bins_b = nb.ang_dist_second() #Angular distribution function upto second neighbor
ddf, bins_d = nb.get_ddf() #Dihedral distribution function upto first neighbor
import matplotlib
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

the_grid = GridSpec(2, 2)
plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(16,14))

plt.subplot(the_grid[0, 0])
plt.title('(a) RDF')
plt.plot(bins_rdf, rdf)
plt.xlabel(r'Distance bins ($\AA$)')

plt.subplot(the_grid[0, 1])
plt.title('(b) ADF-a')
plt.plot(bins_a[:-1], adfa)
plt.xlabel(r'Angle bins ($^\circ$)')

plt.subplot(the_grid[1, 0])
plt.title('(c) ADF-b')
plt.plot(bins_b[:-1], adfb)
plt.xlabel(r'Angle bins ($^\circ$)')

plt.subplot(the_grid[1, 1])
plt.title('(d) DDF')
plt.plot(bins_d[:-1], ddf)
plt.xlabel(r'Angle bins ($^\circ$)')
plt.tight_layout()
```

### How to get XRD paterns

X-ray diffraction patterns act as one of the most important experimental
methods for determining atomic structure. Using Cu-K alpha wavelength,
the theoretical XRD patterns (two-theta and d_hkl dependence) for Si
class above can be obatined as follows.

``` python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
a, b, c = XRD().simulate(atoms=atoms)

the_grid = GridSpec(1,2)
plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(10,5))

plt.subplot(the_grid[0])
plt.bar(a,c)
plt.xlabel('2$\Theta$')
plt.ylabel('XRD intensity')
plt.subplot(the_grid[1])
plt.bar(a,b)
plt.xlabel('d$_{hkl}$')
plt.ylabel('XRD intensity')
plt.tight_layout()
```

### How to make defects

While the above Si atomic structure generated above is perfect/defect
free, in reality there can be several defects present in an atomic
structure such as point defects (vacancies, interstitials,
substituions), line defects (dislocations), surface-defects
(free-surfaces, grain boundaries, stacking faults, interfaces),
volume-defects (voids/pores) etc.

An example of creating vacancy structures using unique Wycoff positions
is shown below:

``` python
from jarvis.analysis.defects.vacancy import Vacancy
#enforces cell-size to be close to 10 Angstroms
vacs = Vacancy(atoms=Si).generate_defects(enforce_c_size=10.0)
len(vacs), vacs[0].to_dict()["defect_structure"].num_atoms
# We find that there are only one unique point vacanc available based on Wycoff-position information
```

Similarly, an example of creating, free surfaces is shown below:

``` python
from jarvis.analysis.defects.surface import wulff_normals, Surface

# Let's create (1,1,1) surface with three layers, and vacuum=18.0 Angstrom
# We center it around origin so that it looks good during visualization
surface_111 = (
    Surface(atoms=Si, indices=[1, 1, 1], layers=3, vacuum=18)
        .make_surface()
        .center_around_origin()
)
print(surface_111)
```

While the above example makes only one surface (111), we can ask
jarvis-tools to provide all symmetrically distinct surfaces as follows:

``` python
from jarvis.analysis.structure.spacegroup import (
    Spacegroup3D,
    symmetrically_distinct_miller_indices,
)
spg = Spacegroup3D(atoms=Si)
cvn = spg.conventional_standard_structure
mills = symmetrically_distinct_miller_indices(max_index=3, cvn_atoms=cvn)
for i in mills:
    surf = Surface(atoms=Si, indices=i, layers=3, vacuum=18).make_surface()
    print ('Index:', i)
    print (surf)
```

Heterostructures of a film and a substrate can be created using ZSL
algorithm as shown in the following example:

``` python
from jarvis.analysis.interface.zur import ZSLGenerator, mismatch_strts, get_hetero, make_interface
film = Surface(atoms=Si, indices=[1, 1, 1], layers = 3, vacuum = 18 ).make_surface().center_around_origin() 
substrate = Surface(atoms=Si, indices=[1, 1, 1], layers = 3, vacuum = 18 ).make_surface().center_around_origin()  
info = make_interface(film=film, subs=substrate)['interface'].center(vacuum=18)
print (info)
```

## How to setup/analyze DFT calculations using VASP

The Vienna Ab initio Simulation Package (VASP) is a package for
performing ab initio quantum mechanical calculations using either
Vanderbilt pseudopotentials, or the projector augmented wave method, and
a plane wave basis set. Manual for VASP is available at:
<https://www.vasp.at/wiki/index.php/The_VASP_Manual> .

Running a VASP calculation requires the following files: `INCAR`,
`POSCAR`, `KPOINTS`, `POTCAR` as well as additional files such as
`vdw_kernel.bindat` for specific types of calculations. While setting up
calculations for one or a few systems/setups should be straight forward,
setting up calculations for thousands of materials and most importantly
making a database out of all those calculations require automated
calculations script collections such as JARVIS-Tools.

Gievn an atomic structure in 1) `jarvis.core.Atoms` format, JARVIS-Tools
2) prepares input files such as `INCAR` etc. as mentioned above and 3)
submits the calculations to your queuing system such as SLURM/PBS using
`jarvis.tasks.vasp` and `jarvis.tasks.queue_jobs`. After a calculations
get completed, 4) automated analysis can be carried out and plots and
webpages are generated. The input file generation and output file
parsing modules for VASP can be found in `jarvis.io.vasp.inputs` and
`jarvis.io.vasp.outputs` modules. The automated analyis and XML
generation for webpages can be found in `jarvis.db.vasp_to_xml` module.
After the xml page creation they are converted using html using XSLT
scripts.

Additionally, a JSON file is created with metadata from all the XML
pages for thousands of materials to easily use in data-analytics/machine
learning applications.The JARVIS-DFT
(<https://jarvis.nist.gov/jarvisdft/>) database primarily uses such a
workflow. Make sure `VASP_PSP_DIR` is declared as a PATH to VASP
pseudopotential directory i.e.

``` bash
$ export VASP_PSP_DIR=YOUR_PATH_TO_PSUEDOPTENTIALS
```

in your ~/.bashrc file.

### How to setup a single calculation

We start by setting up and submitting a single VaspJob:

``` python
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
```

The job.py can now be run on a cluster or on a PC as a python script.
For running this job on a PBS cluster,

``` python
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
```

### How to setup high-throughput calculations

Currently, JARVIS-Tools can be used to submit job with SLURM and PBS
clusters only. For high-throughput automated submissions one can use
pre-build `JobFactory` module that allows automatic calculations for a
series of properties.

``` python
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
```

We provide modules to convert the calculation informato to `XML` which
can be converted to `HTML` using `XSLT`. An example is give below:

``` python
from jarvis.db.vasp_to_xml import VaspToApiXmlSchema
from jarvis.db.restapi import Api
folder="jarvis/jarvis/examples/vasp/SiOptB88vdW"
filename = "JVASP-1002.xml"
VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)
```

### How to plot electronic bandstructure and DOS

If you use the workflow used above, the density of states plot can be
obtained using thr `vasprun.xml` file in MAIN-RELAX folder while the
band-structure plot is obtained using `vasprun.xml` in MAIN-BAND folder.

``` python
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
```

### How to obtain elastic constants

### How to plot generate an STM/STEM image

### How to plot generate a dielectric function spectra and solar eff.

### How to generate/use electronic Wannier tight binding model

### How to generate Fermi-surfaces

### How to run BoltzTrap for transport properties

### How to make heterostructures/interfaces

### How to get IR/Raman spectra

### How to get piezoelectic/dielecrric/BEC constants

### How to get electric field gradients

### How to get work-function of a surface

### How to get exfoliation energy of a 2D material

## How to run/analyze MD static/dynamic calculation using LAMMPS

Molecular dynamics/classical force-field calculations can be carried out
with LAMMPS software as in JARVIS-FF. An example for running LAMMPS is
given below. Here, a `LammpsJob` module is defined with the help of
atoms, pair-style, coefficient, and template file (\*.mod file) to
control the calculations.

### How to run calculation

``` python
from jarvis.tasks.lammps.lammps import LammpsJob, JobFactory
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from jarvis.analysis.structure.spacegroup import Spacegroup3D


# atoms = Atoms.from_poscar('POSCAR')
# Get Aluminum FCC from JARVIS-DFT database
tmp_dict = get_jid_data(jid="JVASP-816", dataset="dft_3d")["atoms"]
atoms = Atoms.from_dict(tmp_dict)

# Get conventional cell
spg = Spacegroup3D(atoms)
cvn_atoms = spg.conventional_standard_structure

# Set-up path to force-field/potential file, .mod file. and lammps executable
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

# Test in a high-throughput
job_fact = JobFactory(pair_style="eam/alloy", name="my_first_lammps_run")
job_fact.all_props_eam_alloy(atoms=cvn_atoms, ff_path=ff, lammps_cmd=cmd)
```

### How to analyze data

An example to parse LAMMPS calculation folder using the above workflow
is shown below:

``` python
from jarvis.io.lammps.outputs import parse_material_calculation_folder
folder = '/home/users/knc6/Software/jarvis/jarvis/examples/lammps/Aleam'
data = parse_material_calculation_folder(folder)
print (data)
```

The calculation data can now be converted into XML files as follows. The
XML with the help of XSLT is converted into an HTML page.

``` python
from jarvis.db.lammps_to_xml import write_xml
write_xml(data=data,filename='JLMP-123.xml')
```

## How to run/analyze DFT static calculation using Quantum espresso

Quantum ESPRESSO is a suite for first-principles electronic-structure
calculations and materials modeling, distributed for free and as free
software under the GNU General Public License. It is based on
density-functional theory, plane wave basis sets, and pseudopotentials.

### How to setup a single calculation

An example for running QE simulation is shown below:

``` python
from jarvis.core.kpoints import Kpoints3D
from jarvis.core.atoms import Atoms
box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
elements = ["Si", "Si"]
Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
print(Si)
kp = Kpoints3D().automatic_length_mesh(
    lattice_mat=Si.lattice_mat, length=20
)
qe = QEinfile(Si, kp)
qe.write_file()
kp = Kpoints3D().kpath(atoms=Si)
qe = QEinfile(Si, kp)
qe.write_file("qe.in2")
sp = qe.atomic_species_string()
sp = qe.atomic_cell_params()
print("sp", sp)
print(qe.input_params['system_params']['nat'])
$PATH_TO_PWSCF/pw.x -i qe.in
```

### How to setup high-throughput calculations

## How to traing JARVIS-CFID ML models using sklearn/lightgbm

There are several methods to train atomistic property ML models such as
based on hand-crafted descritprs and graph neural network. Examples of
such methods are: JARVIS-CFID (Classical Force-Field Inspired
Descriptors) for descriptors based training and JARVIS-ALIGNN (Atomistic
Line Graph Neural Network) based on GNNs. In this section we discuss the
JARVIS-CFID ( `jarvis.ai.descriptors.cfid`), which can be used for
training models with only chemical formula or chemical formula+structure
information.

### How to train chemical formula only datasets

For each chemical formula, we can obtain <span class="title-ref">438
descriptors</span> consisting of features such as avergae
electronegativity, average boiling points of elements etc. An example of
getting descriptors isshown below:

``` python
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie
from jarvis.ai.pkgs.lgbm.regression import regression
from jarvis.ai.descriptors.cfid import get_chem_only_descriptors

# Load a dataset, you can use pandas read_csv also to generte my_data
# Here is a sample dataset
my_data = [
    ["CoAl", 1],
    ["CoNi", 2],
    ["CoNb2Ni5", 3],
    ["Co1.2Al2.3NiRe2", 4],
    ["Co", 5],
    ["CoAlTi", 1],
    ["CoNiTi", 2],
    ["CoNb2Ni5Ti", 3],
    ["Co1.2Al2.3NiRe2Ti", 4],
    ["CoTi", 5],
    ["CoAlFe", 1],
    ["CoNiFe", 2],
    ["CoNb2Ni5Fe", 3],
    ["Co1.2Al2.3NiRe2Fe", 4],
    ["CoFe", 5],
]


# Convert my_data to numpy array
X = []
Y = []
IDs = []
for ii, i in enumerate(my_data):
    X.append(get_chem_only_descriptors(i[0]))
    Y.append(i[1])
    IDs.append(ii)

X = np.array(X)
Y = np.array(Y).reshape(-1, 1)
IDs = np.array(IDs)
```

Now, we can use different ML algorithms on the descriptors and dataset
such as linear regression, random forest, gradient boosting etc.

An example, for using LightGBM with jarvis-tools wrapper code is shown
below:

``` python
# Train a LightGBM regression model
config = {"n_estimators": 5, "learning_rate": 0.01, "num_leaves": 2}
# The regression module does feature pre-processing as well
# Change config settings to improve model such as by hyper-parameter tuning
info = regression(X=X, Y=Y, jid=IDs, feature_importance=False, config=config)


# Print performance metrices
# Print performance metrices
print(
    'r2=',info["reg_scores"]["r2"],
    'MAE=',info["reg_scores"]["mae"],
    'RMSE=',info["reg_scores"]["rmse"],
)
```

### How to train regression model

Suppose we have 60000 materials, and we get 1557 descriptor for each
material (438 chemical as above as well as structure and charge
descriptors), we will have a 60000x1557 matrix. Let's call this matrix
as 'x' or input matrix. Next, we can get target ('y') data either from
DFT, FF calculations or experiments. For example, we can choose
formation energies of 60000 materials in the JARVIS-DFT as the dtarget
data giving 60000x1 matrix.

Now, we can use a ML/AI algorithm to establish statistical relation
between the x and y data. Once trained we get a trained model, which can
be stored in say pickle or joblib format.

For a new material now, it can be converted into CFID i.e. 1x1557 matrix
which when fed to the model will give 1x1 prediction hence the ML
prediction. We can use a range of ML algorithms such as linear
regression, decision trees, Gaussian processes etc. We find with CFID
descriptors, gradient boosting decision trees (especially in LightGBM)
gives one of the most accurate results. We provide tools to run with
major ML packages such as scikit-learn, tensorflow, pytorch, lightgbm
etc. Example-1:

``` python
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

## How to traing JARVIS-ALIGNN ML models using PyTorch

### How to train regression model

How to train classification model ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

## How to use quantum computation algorithms using Qiskit/Tequila/Pennylane

Quantum chemistry is one of the most attractive applications for quantum
computations. Predicting the energy levels of a Hamiltonian is a key
problem in quantum chemistry. Variational quantum eigen solver (VQE) is
one of the most celebrated methods for predicting an approximate ground
state of a Hamiltonian on a quantum computer following the variational
principles of quantum mechanics.VQE utilizes Ritz variational principle
where a quantum computer is used to prepare a wave function ansatz of
the system and estimate the expectation value of its electronic
Hamiltonian while a classical optimizer is used to adjust the quantum
circuit parameters in order to find the ground state energy. A typical
VQE task is carried out as follows: an ansatz/circuit model with tunable
parameters is constructed and a quantum circuit capable of representing
this ansatz is designed. In this section, we show a few examples to
apply quantum algorithms for solids using Wannier-tight binding
Hamiltonians (WTBH). WTBHs can be generated from several DFT codes.
Here, we use JARVIS-WTBH database.

### How to generate circuit model

Developing a heuristic quantum circuit model is probably the most
challenging part of applying quantum algorithms. Fortunately, there are
few well-known generalized models that we can use or generate ourselves.
There are several circuit models (for a fixed number of qubits and
repeat units) available in `jarvis.core.circuits.`. In the following
example, we use circuit6/EfficientSU2 model and use it to predict
electronic energy levels (at a K-point in the Brillouin zone) of FCC
Aluminum using a WTBH.

``` python
from jarvis.db.figshare import get_wann_electron, get_wann_phonon, get_hk_tb
from jarvis.io.qiskit.inputs import HermitianSolver
from jarvis.core.circuits import QuantumCircuitLibrary
from qiskit import Aer

backend = Aer.get_backend("statevector_simulator")
# Aluminum JARVIS-ID: JVASP-816
wtbh, Ef, atoms = get_wann_electron("JVASP-816") 
kpt = [0.5, 0., 0.5] # X-point
hk = get_hk_tb(w=wtbh, k=kpt)
HS = HermitianSolver(hk)
n_qubits = HS.n_qubits()
circ = QuantumCircuitLibrary(n_qubits=n_qubits).circuit6()
en, vqe_result, vqe = HS.run_vqe(var_form=circ, backend=backend)
vals,vecs = HS.run_numpy()
# Ef: Fermi-level
print('Classical, VQE (eV):', vals[0]-Ef, en-Ef)
print('Show model\n', circ)
```

### How to run cals. on simulators

In the above example, we run simulations on `statevector_simulator`.
Qiskit provides several other simulators, which can also be used.

### How to run cals. on actual quantum computers

To run calculations on real quantum computers, we just replace the
`backend` parameter above such as the following:

``` python
token='Get Token from your IBM account' 
qiskit.IBMQ.save_account(token)
provider = IBMQ.load_account()
backend = provider.get_backend('ibmq_5_yorktown')
```

Your job will put in a queue and as the simulation complete result will
be sent back to you. Note that there might be a lot of jobs in the queue
already, so it might take a while. You may run simulations using IBM GUI
or use something like Jupyter notebook/Colab notebook.
