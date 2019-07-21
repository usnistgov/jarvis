
from __future__ import unicode_literals, print_function

"""
Helper function for running LAMMPS
Used for defects, surface and phonon calculations
"""

from monty.json import MontyEncoder
from numpy import matrix
import shutil, time
from ase import *
import numpy as np

try:
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS, write_FORCE_CONSTANTS
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
except:
    pass
import glob, fileinput
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.core.surface import (
    Slab,
    SlabGenerator,
    generate_all_slabs,
    get_symmetrically_distinct_miller_indices,
)
from ase.lattice.surface import surface
from pymatgen.ext.matproj import MPRester
import operator
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from jarvis.point_defect.vacancy import vac_antisite_def_struct_gen, chempot_struct
from jarvis.plane_defect.surface import pmg_surfer, surfer
import numpy as np, time, json
import sys, os, subprocess, socket
from pymatgen.io.ase import AseAtomsAdaptor
from jarvis.tools.zipit import ZipDir
from jarvis.phonopy.phonon import get_phonopy_atoms
import fortranformat as fform
from pymatgen.core.structure import Structure
from ase.calculators.lammpsrun import LAMMPS, prism
try:
    input_nobox = os.environ["input_nobox"]
    input_box = os.environ["input_box"]
    lammps_exec = os.environ["lammps_exec"]
except:
    pass


def write_lammps_data(structure=None, file="", write_tmp_file=True):
    """
        write lammps structure data
        from ase with custom modifications

        Args:
            structure: Structure object
            file:  intended file to write in
        """
    structure.sort()
    if write_tmp_file == True:
        structure.to(fmt="poscar", filename="new_pymatgen_slab.vasp")
    atoms = AseAtomsAdaptor().get_atoms(structure)

    f = open(file, "w")
    f.write("datafile (written by JARVIS-FF) \n\n")
    symbols = atoms.get_chemical_symbols()
    import ase.io.vasp

    # ase.io.vasp.write_vasp("POSCAR.1x1x1",atoms, label='444supercell',direct=True,sort=False)
    n_atoms = len(symbols)
    f.write("%d \t atoms \n" % n_atoms)
    species = [tos for tos in Poscar(structure).site_symbols]
    # species = [tos.symbol
    #               for tos in structure.types_of_specie]
    n_atom_types = len(species)
    print("species", species)
    f.write("%d  atom types\n" % n_atom_types)
    p = prism(atoms.get_cell())
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()
    f.write("0.0 %s  xlo xhi\n" % xhi)
    f.write("0.0 %s  ylo yhi\n" % yhi)
    f.write("0.0 %s  zlo zhi\n" % zhi)
    f.write("%s %s %s  xy xz yz\n" % (xy, xz, yz))
    f.write("\n\n")
    f.write("Atoms \n\n")
    for i, r in enumerate(map(p.pos_to_lammps_str, atoms.get_positions())):
        c = 0.0
        s = species.index(symbols[i]) + 1

        f.write("%6d %3d %6f %s %s %s\n" % ((i + 1, s, c) + tuple(r)))
    f.close()


def write_lammps_in(
    structure=None,
    lammps_in="init.mod",
    lammps_in1="potential.mod",
    lammps_in2="in.main",
    lammps_trj=None,
    lammps_data=None,
    parameters={},
):
    """
        write lammps input file
        from ase with custom modifications
        LAMMPS input is devided into three parts

        Args:
            structure: Structure object
            lammps_in: generally"init.mod", with unit and conversion factor information
            lammps_in1: generally "potential.mod", with force-field/potential style and element tyoe information
            lammps_in2: generally "in.elastic", a generic main input file to be fed in LAMMPS usin lmp_*<...,parameters['exec']
            parameters: input parameters
            
        """
    structure.sort()
    f = open(lammps_in, "w")
    f1 = open(lammps_in1, "w")  # potential.mod
    f2 = open(lammps_in2, "w")
    f.write(
        ('variable dump_file string "%s"\n' % lammps_trj)
        + ("variable up  equal 1.0e-6\n")
        + ("variable cfac  equal 1.0e-4\n")
        + ("variable cunits  string GPa\n")
        + ("variable etol  equal 0\n")
        + ("variable ftol  equal 1.0e-10\n")
        + ("variable maxiter equal 1000\n")
        + ("variable maxeval equal 10000\n")
        + ("variable dmax equal 1.0e-2\n")
        + ('variable data_file string "%s"\n' % "data")
    )
    atoms = AseAtomsAdaptor().get_atoms(structure)
    pbc = atoms.get_pbc()
    if "control_file" in parameters:
        f2.write("include %s \n" % parameters["control_file"])
    if "units" in parameters:
        f.write("units %s \n" % parameters["units"])
    else:
        f.write("units metal \n")
    if "atom_style" in parameters:
        f.write("atom_style %s \n" % parameters["atom_style"])
    else:
        f.write("atom_style atomic \n")
    if "boundary" in parameters:
        f.write("boundary %s \n" % parameters["boundary"])
    else:
        f.write("boundary %c %c %c \n" % tuple("sp"[x] for x in pbc))
    f.write("atom_modify sort 0 0.0 \n")
    for key in ("neighbor", "newton"):
        if key in parameters:
            f.write("%s %s \n" % (key, parameters[key]))
    f.write("\n")
    # If no_data_file,
    # write the simulation box and the atoms
    species = [tos for tos in Poscar(structure).site_symbols]
    # species = [tos.symbol for tos in structure.types_of_specie]
    n_atom_types = len(species)
    species_i = dict([(s, i + 1) for i, s in enumerate(species)])
    f.write("read_data %s\n" % "data")
    # interaction
    f.write("\n### interactions \n")
    if "lib" in parameters:
        lib = parameters["lib"]
        f1.write("%s \n" % lib)
    if ("pair_style" in parameters) and ("pair_coeff" in parameters):
        pair_style = parameters["pair_style"]
        f1.write("pair_style %s \n" % pair_style)
        symbols = atoms.get_chemical_symbols()
        # species = [tos.symbol.replace("Mo","M") for tos in structure.types_of_specie] #For REBO Mo-S
        # species = [tos.symbol for tos in structure.types_of_specie]
        print("site symbolss", Poscar(structure).site_symbols)
        species = [tos for tos in Poscar(structure).site_symbols]
        if parameters["pair_style"] == "rebomos":

            species = [tos.replace("Mo", "M") for tos in Poscar(structure).site_symbols]
        tag = ""
        for i in species:
            tag = tag + " " + i
        pair_coef = "* * " + str(parameters["pair_coeff"]) + " " + tag
        f1.write("pair_coeff %s \n" % pair_coef)

        masses = []
        for m in atoms.get_masses():
            if m not in masses:
                masses.append(m)
        count = 0
        for i in masses:
            count = count + 1
            f.write("mass" + " " + str(count) + " " + str(i) + "\n")
    else:
        # default interaction
        f.write("pair_style lj/cut 2.5 \n" + "pair_coeff * * 1 1 \n" + "mass * 1.0 \n")
    f1.write("neighbor 1.0 nsq\n")
    f1.write("neigh_modify once no every 1 delay 0 check yes\n")
    if "min" not in parameters:
        f1.write("min_style  cg\n")
        f1.write("min_modify           dmax ${dmax} line quadratic\n")
    f1.write("thermo          1\n")
    f1.write(
        "thermo_style custom step temp press cpu pxx pyy pzz pxy pxz pyz ke pe etotal vol lx ly lz atoms\n"
    )
    # f1.write('thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol\n' )
    f1.write("thermo_modify norm no\n")
    #   if 'thermo_style' in parameters:
    #       f.write('thermo_style %s\n' % parameters['thermo_style'])
    #   else:
    #       f.write(('thermo_style custom %s\n') %
    #               (' '.join(self._custom_thermo_args)))
    #   if 'thermo_modify' in parameters:
    #       f.write('thermo_modify %s\n' % parameters['thermo_modify'])
    #   else:
    #       f.write('thermo_modify flush yes\n')
    #   if 'thermo' in parameters:
    #       f.write('thermo %s\n' % parameters['thermo'])
    #   else:
    #       f.write('thermo 1\n')
    #   if 'minimize' in parameters:
    #       f.write('minimize %s\n' % parameters['minimize'])
    #   if 'run' in parameters:
    #       f.write('run %s\n' % parameters['run'])
    #   if not (('minimize' in parameters) or ('run' in parameters)):
    #       f.write('run 0\n')
    #   if 'dump' in parameters:
    #       f.write('dump %s\n' % parameters['dump'])
    #   else:
    #       f.write('dump dump_all all custom 1 %s id type x y z vx vy vz fx fy fz\n' % lammps_trj)
    #   f.write('print __end_of_ase_invoked_calculation__\n')
    #   f.write('log /dev/stdout\n')
    if "fix" in parameters:
        if parameters["fix"]:
            for i in parameters["fix"]:
                f1.write("fix %s\n" % i)
    f.close()
    f1.close()
    f2.close()


def analyz_loge(log="log.lammps"):
    import sys

    """ 
    Analyzes log.lammps file,
    Please note, the output format heavily depends on the input file
    A generic inpu is taken here

    Args:
        log: path to log.lammps file
    Returns:
          en: energy/atom
          press: pressure
          toten: total energy
          cij: elastic constants
    """
    en = 0
    press = 0
    c11 = 0
    c22 = 0
    c33 = 0
    c44 = 0
    c55 = 0
    c66 = 0
    c12 = 0
    c13 = 0
    c23 = 0
    c14 = 0
    c15 = 0
    c16 = 0
    c14 = 0
    c24 = 0
    c25 = 0
    c26 = 0
    c34 = 0
    c35 = 0
    c36 = 0
    c45 = 0
    c46 = 0
    c56 = 0
    try:
        logfile = open(log, "r")
        lines = logfile.read().splitlines()
        for i, line in enumerate(lines):
            if "Loop time of" in line:
                toten = float(lines[i - 1].split()[12])
                press = float(lines[i - 1].split()[2])
                press = float(press) * 0.0001
                en = float(lines[i - 1].split()[12]) / float(lines[i - 1].split()[17])
                break
        logfile.close()
    except:
        pass
    try:
        logfile = open(log, "r")
        lines = logfile.read().splitlines()
        for i, line in enumerate(lines):
            if 'print "Elastic Constant C11all = ${C11all} ${cunits}"' in line:
                c11 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C22all = ${C22all} ${cunits}"' in line:
                c22 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C33all = ${C33all} ${cunits}"' in line:
                c33 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C12all = ${C12all} ${cunits}"' in line:
                c12 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C13all = ${C13all} ${cunits}"' in line:
                c13 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C23all = ${C23all} ${cunits}"' in line:
                c23 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C44all = ${C44all} ${cunits}"' in line:
                c44 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C55all = ${C55all} ${cunits}"' in line:
                c55 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C66all = ${C66all} ${cunits}"' in line:
                c66 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C14all = ${C14all} ${cunits}"' in line:
                c14 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C16all = ${C16all} ${cunits}"' in line:
                c16 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C24all = ${C24all} ${cunits}"' in line:
                c24 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C25all = ${C25all} ${cunits}"' in line:
                c25 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C26all = ${C26all} ${cunits}"' in line:
                c26 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C34all = ${C34all} ${cunits}"' in line:
                c34 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C35all = ${C35all} ${cunits}"' in line:
                c35 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C36all = ${C36all} ${cunits}"' in line:
                c36 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C45all = ${C45all} ${cunits}"' in line:
                c45 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C46all = ${C46all} ${cunits}"' in line:
                c46 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C56all = ${C56all} ${cunits}"' in line:
                c56 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
        logfile.close()
    except:
        pass
    return (
        round(en, 2),
        round(press, 2),
        float(toten),
        round(float(c11), 1),
        round(float(c22), 1),
        round(float(c33), 1),
        round(float(c12), 1),
        round(float(c13), 1),
        round(float(c23), 1),
        round(float(c44), 1),
        round(float(c55), 1),
        round(float(c66), 1),
        round(float(c14), 1),
        round(float(c16), 1),
        round(float(c24), 1),
        round(float(c25), 1),
        round(float(c26), 1),
        round(float(c34), 1),
        round(float(c35), 1),
        round(float(c36), 1),
        round(float(c45), 1),
        round(float(c46), 1),
        round(float(c56), 1),
    )


def read_data(data=None, ff=None):
    """
    Read LAMMPS data file

    Args:
        data: data file path
        ff: potential.mod/potential information file path
    Returns:
          struct: Structure object
    """
    pot_file = open(ff, "r")
    lines = pot_file.read().splitlines()
    symb = []
    count = 0
    from pymatgen.core.periodic_table import Element

    for i, line in enumerate(lines):
        if "pair_coeff" in line.split():
            sp = line.split()
            print("spsplit", sp, os.getcwd())
            for el in sp:
                try:
                    if Element(el):
                        # if el=='M':
                        #    el='Mo'
                        # count=count+1
                        # if count >4:
                        symb.append(el)
                except:
                    pass
    print("symb=", symb)

    f = open(data, "r")
    lines = f.read().splitlines()
    for i, line in enumerate(lines):
        if "atoms" in line.split():
            natoms = int(line.split()[0])
        if "types" in line.split():
            print(line)
            ntypes = int(line.split()[0])
        if "xlo" in line.split():
            xlo = float(line.split()[0])
            xhi = float(line.split()[1])
        if "ylo" in line.split():
            ylo = float(line.split()[0])
            yhi = float(line.split()[1])
        if "zlo" in line.split():
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])
        if "xy" in line.split():
            xy = float(line.split()[0])
            xz = float(line.split()[1])
            yz = float(line.split()[2])
    if len(symb) != ntypes:
        print("Something wrong in atom type assignment", len(symb), ntypes)
        sys.exit()
    lat = Lattice([[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]])
    typ = np.empty((natoms), dtype="S20")
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    z = np.zeros((natoms))
    q = np.zeros((natoms))
    coords = list()
    for i, line in enumerate(lines):
        if "Atoms" in line.split():
            for j in range(0, natoms):
                # print int(((lines[j+2]).split()[1]))-1
                typ[j] = symb[int(((lines[i + j + 2]).split()[1])) - 1]
                q[j] = (lines[i + j + 2]).split()[2]
                x[j] = (lines[i + j + 2]).split()[3]
                y[j] = (lines[i + j + 2]).split()[4]
                z[j] = (lines[i + j + 2]).split()[5]
                coords.append([x[j], y[j], z[j]])
    f.close()
    # print ("info",(typ),'coo',(coords),'latt',lat)
    pot_file.close()
    typ_sp = [str(i, "utf-8") for i in typ]
    struct = Structure(lat, typ_sp, coords, coords_are_cartesian=True)
    # print struct
    # finder = SpacegroupAnalyzer(struct)
    # num=finder.get_spacegroup_symbol()
    # print(num)
    return struct


def smart_vac(strt=None, parameters=None):
    """
    Function to get all vacancy formation energies

    Args:
        strt: Structure object
        parameters: parameters with LAMMPS inputs
    """
    parameters["control_file"] = input_nobox  #'/users/knc6/inelast_nobox.mod'
    vac_arr = []
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()

    cellmax = 2  # int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    ase_atoms = ase_atoms * (cellmax, cellmax, cellmax)
    # if len(ase_atoms) >200:
    #   cellmax=1
    # else:
    #   cellmax=2
    # cellmax=int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    print("type of trt is= celmmax", type(strt), cellmax)
    try:
        print("int(strt.composition.num_atoms)", int(strt.composition.num_atoms))
        print(int(strt.ntypesp))
    except:
        pass
    # cellmax=int(strt.composition.num_atoms)+int(strt.ntypesp)
    vac_done = 0
    tol = 0.01  # change 0.1
    try:
        vac = vac_antisite_def_struct_gen(cellmax=cellmax, struct=strt)
    except:
        print("Failed for v_1", os.getcwd())
        pass
    def_list = [100000 for y in range(len(vac) - 1)]
    while vac_done != 1:
        try:
            vac = vac_antisite_def_struct_gen(cellmax=cellmax, struct=strt)
        except:
            print("Failed for v_2", os.getcwd())
            pass
        if vac not in vac_arr:
            try:
                vac_arr.append(vac)
                print(
                    "in smart_vac(strt=None), cellmax,vac,vac_done=",
                    cellmax,
                    vac,
                    vac_done,
                )
                def_list2, header_list = def_energy(vac=vac, parameters=parameters)
                diff = matrix(def_list) - matrix(def_list2)
                diff_arr = np.array(diff).flatten()
                print("in smart_vac(strt=None diff_arr=", diff_arr)
                print("sleepig for 5")
                # time.sleep(5)
            except:
                print("Failed for v_3", os.getcwd())
                pass
            if any(diff_arr) > tol:
                # for el in diff_arr:
                #     if abs(el)>tol :
                # print ("in smart_vac(strt=None abs_el=",abs(el))
                vac_done = 0
                cellmax = cellmax + 1
                ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
                ase_atoms = ase_atoms * (cellmax, cellmax, cellmax)
                if len(ase_atoms) > 50:
                    vac_done = 1
                    break
                def_list = def_list2
            else:
                vac_done = 1
    #        cellmax=cellmax+1
    return def_list, header_list


def smart_surf(strt=None, parameters=None, layers=3, tol=0.5):
    """
    Function to get all surface energies

    Args:
        strt: Structure object
        parameters: parameters with LAMMPS inputs
        layers: starting number of layers
        tol: surface energy tolerance for convergence
    Returns:
          surf_list: list of surface energies
          surf_header_list: list of surface names

    """
    parameters["control_file"] = input_nobox  #'/users/knc6/inelast_nobox.mod'
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    layers = 3
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    for i in indices:
        ase_slab = surface(ase_atoms, i, layers)
        ase_slab.center(vacuum=15, axis=2)
        if len(ase_slab) < 50:
            layers = 3
    surf_arr = []
    surf_done = True
    try:
        surf = surfer(mat=strt, layers=layers)
        surf_list = [100000 for y in range(len(surf) - 1)]
        print("in smart_surf :surf,surf_list=", surf, surf_list)
    except:
        print("Failed at s1", os.getcwd())
        pass
    while surf_done:
        layers = layers + 1
        indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
        ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
        for i in indices:
            ase_slab = surface(ase_atoms, i, layers)
            ase_slab.center(vacuum=15, axis=2)
            # if len(ase_slab) > 100:
            #   surf_done=True
            # if (ase_slab.get_cell()[2][2]) > 40:
            #   surf_done=True

        try:
            surf = surfer(mat=strt, layers=layers)
        except:

            print("Failed at s2", os.getcwd())
            pass
        if surf not in surf_arr:
            surf_arr.append(surf)
            try:
                surf_list2, surf_header_list = surf_energy(
                    surf=surf, parameters=parameters
                )
                print("in smart_surf :surf2,surf_list2=", surf_list2, surf_header_list)
                diff = matrix(surf_list) - matrix(surf_list2)
                print(
                    "in smart_surf :surf3,surf_list3=",
                    matrix(surf_list),
                    matrix(surf_list2),
                )
                diff_arr = np.array(diff).flatten()
            except:
                print("Failed for s_3", os.getcwd())
                pass
            if len(ase_slab) > 50:
                surf_done = True
                break
                # print ("layersssssssssssssssssssssssss",layers,surf_done)
                break
            if any(diff_arr) > tol:
                # for el in diff_arr:
                #    if abs(el)>tol :
                #        print ("in smart_surf :abs el=",abs(el))
                surf_done = True
                surf_list = surf_list2
            else:
                surf_done = False
    return surf_list, surf_header_list


def surf_energy(surf=[], parameters={}):
    """
    Function to get specific surface energies

    Args:
        surf: list of syrfaces
        parameters: LAMMPS input parameters
    Returns:
          surf_list: list of surface energies
          surf_header_list: list of surface names
    """
    # def surf_energy(surf=[],parameters = {'pair_style':'eam/alloy','pair_coeff':'/scratch/lfs/kamal/POSMAT/Automatic2/Al03.eam.alloy','atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'}):
    surf_list = []
    surf_header_list = []
    fi = str("SURF.INFO")
    f = open(fi, "w")
    for s in surf:
        enp, strt, forces = run_job(mat=s, parameters=parameters)
        # strt=Structure.from_file(contc)
        m = strt.lattice.matrix
        if s.comment.split("@")[0] == "sbulk":

            gs_energy = enp
            print("in surf_energy gs_energy for", s.comment, gs_energy)
        else:
            surf_en = (
                16
                * (
                    float(enp) * (strt.composition.num_atoms)
                    - float(strt.composition.num_atoms) * float(gs_energy)
                )
                / (2 * np.linalg.norm(np.cross(m[0], m[1])))
            )
            print("in surf_energy surf_en for", s.comment, surf_en)
            surf_list.append(surf_en)
            print(s.comment, "=", surf_en)
            line = str(s.comment) + str("=") + str(surf_en) + "\n"
            f.write(line)
    f.close()
    return surf_list, surf_header_list


def get_chem_pot(s1=None, s2=None, parameters={}):
    """
    Get chemical potential given perfect and defect structures

    Args:
        s1: perfect Structure object
        s2: defect Structure object
        parameters: LAMMPS input parameter
    Returns:
           enp: energy per atom
    """
    # def get_chem_pot(s1=None,s2=None,parameters= {'pair_style':'eam/alloy','pair_coeff':'/scratch/lfs/kamal/POSMAT/Automatic2/Al03.eam.alloy','atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'}):
    s1.sort()
    s2.sort()
    s3 = (set(s2)).symmetric_difference(set(s1))  # list(set(s1)-set(s2))
    # s3=[]
    # for el1 in s1:
    #     for el2 in s2:
    #         if el1 !=el2 and el1 not in s3:
    #            s3.append(el1)
    #         if el1 !=el2 and el2 not in s3:
    #            s3.append(el2)
    #
    # from pymatgen.analysis.structure_matcher import StructureMatcher
    # print "Mather",StructureMatcher().fit(s1, s2)
    print("s3   is   ", type(s1), type(s2), s3)
    uniq = []
    for q in s3:
        el = q._species.elements
        for j in el:
            print("j is ", j)
            if j not in uniq:
                uniq.append(j)
                a, b = chempot_struct(j.symbol)
                # a, b = get_struct_from_jv(j)
                p = Poscar(b)
                p.comment = str(a)
                enp, strt, forces = run_job(mat=p, parameters=parameters)
    if len(uniq) > 1:
        print("uniq problem", uniq)
    return enp


def calc_forces(mat=None, parameters={}):
    """
    Calculate forces on atoms
   
    """
    enp, strt, forces = run_job(mat=mat, parameters=parameters)
    return forces


def do_phonons(strt=None, parameters=None, c_size=15):
    """
    Setting up phonopy job using LAMMPS

    Args:
        strt: Structure object
        parameters: LAMMPS input file parameters
        c_size: cell-size 
    
    """
    # spg_strt = SpacegroupAnalyzer(strt).get_conventional_standard_structure()
    p = get_phonopy_atoms(mat=strt)
    bulk = p

    dim1 = int((float(c_size) / float(max(abs(strt.lattice.matrix[0]))))) + 1
    dim2 = int(float(c_size) / float(max(abs(strt.lattice.matrix[1])))) + 1
    dim3 = int(float(c_size) / float(max(abs(strt.lattice.matrix[2])))) + 1
    Poscar(strt).write_file("POSCAR")
    tmp = strt.copy()
    tmp.make_supercell([dim1, dim2, dim3])
    Poscar(tmp).write_file("POSCAR-Super.vasp")

    # phonon = Phonopy(get_phonopy_atoms(tmp), [[dim1, 0, 0], [0, dim2, 0], [0, 0, dim3]])  # ,
    phonon = Phonopy(bulk, [[dim1, 0, 0], [0, dim2, 0], [0, 0, dim3]])  # ,
    print("[Phonopy] Atomic displacements1:", bulk)
    print("[Phonopy] Atomic displacements2:", phonon, dim1, dim2, dim3)
    phonon.generate_displacements(distance=0.03)
    disps = phonon.get_displacements()
    print("[Phonopy] Atomic displacements3:", disps)
    for d in disps:
        print("[Phonopy]", d[0], d[1:])
    supercells = phonon.get_supercells_with_displacements()

    # Force calculations by calculator
    set_of_forces = []
    disp = 0
    for scell in supercells:
        cell = Atoms(
            symbols=scell.get_chemical_symbols(),
            scaled_positions=scell.get_scaled_positions(),
            cell=scell.get_cell(),
            pbc=True,
        )
        disp = disp + 1

    mat = Poscar(AseAtomsAdaptor().get_structure(cell))
    mat.comment = str("disp-") + str(disp)
    parameters["min"] = "skip"
    parameters["control_file"] = "/users/knc6/in.phonon"
    # a,b,forces=run_job(mat=mat,parameters={'min':'skip','pair_coeff': '/data/knc6/JARVIS-FF-NEW/ALLOY4/Mishin-Ni-Al-2009.eam.alloy', 'control_file': '/users/knc6/in.phonon', 'pair_style': 'eam/alloy', 'atom_style': 'charge'})
    a, b, forces = run_job(mat=mat, parameters=parameters)
    print("forces=", forces)
    drift_force = forces.sum(axis=0)
    print("drift forces=", drift_force)
    print("[Phonopy] Drift force:", "%11.5f" * 3 % tuple(drift_force))
    # Simple translational invariance
    for force in forces:
        force -= drift_force / forces.shape[0]
    set_of_forces.append(forces)
    phonon.produce_force_constants(forces=set_of_forces)

    write_FORCE_CONSTANTS(phonon.get_force_constants(), filename="FORCE_CONSTANTS")
    print()
    print("[Phonopy] Phonon frequencies at Gamma:")


# for i, freq in enumerate(phonon.get_frequencies((0, 0, 0))):
#    print ("[Phonopy] %3d: %10.5f THz" %  (i + 1, freq) # THz)


def def_energy(vac=[], parameters={}):
    """
    Get specific defect formation energy

    Args:
         vac: vacancy structures
         parameters: LAMMPS input parameters
    Returns:
           def_list: defect energy list
           header_list: defect names
    """
    # def def_energy(vac=[],parameters={'pair_style':'eam/alloy','pair_coeff':'/scratch/lfs/kamal/POSMAT/Automatic2/Al03.eam.alloy','atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'}):
    def_list = []
    header_list = []
    fi = str("DEF.INFO")
    f = open(fi, "w")
    for v in vac:
        enp, strt, forces = run_job(mat=v, parameters=parameters)
        # enp,contc=run_job(mat=v,incar=incar,kpoints=kpoints,jobname=str(v.comment))
        # strt=Structure.from_file(contc)
        comm = str(v.comment)
        print("running def_energy for =", comm)
        header_list.append(comm)
        if comm.split("@")[0] == "bulk":
            gs_energy = float(enp) * float(strt.composition.num_atoms)
            gs_str = strt
            print("in def_energy gs_energy for", comm, gs_energy)
        else:
            print("strt", strt)
            print("gs_str", gs_str)
            chem_pot = get_chem_pot(gs_str, strt, parameters=parameters)
            if comm.startswith("intl"):
                chem_pot = 0.0 - float(chem_pot)
            def_en = (
                (float(enp) * float(strt.composition.num_atoms))
                - float(gs_energy)
                + chem_pot
            )
            print("in def_energy def_en for", comm, def_en, "chem_pot", chem_pot)
            def_list.append(def_en)
            print(v.comment, "=", def_en)
            line = str(v.comment) + str("=") + str(def_en) + "\n"
            f.write(line)
    f.close()
    return def_list, header_list


# def main(p=None, parameters={'pair_style':'rebomos','pair_coeff':'/scratch/lfs/kamal/JARVIS/All2/MoS2/MoS.REBO.set5b','atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'}):
def main(p=None, parameters={}, c_size=10):
    """
    Master function to run LAMMPS job

    Args:
        p: Poscar object
        parameters: LAMMPS input parameters
        c_size:cell size
    """
    # p=Poscar.from_file("POSCAR")
    sg_mat = SpacegroupAnalyzer(p.structure)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    dim1 = int((float(c_size) / float(max(abs(mat_cvn.lattice.matrix[0]))))) + 1
    dim2 = int(float(c_size) / float(max(abs(mat_cvn.lattice.matrix[1])))) + 1
    dim3 = int(float(c_size) / float(max(abs(mat_cvn.lattice.matrix[2])))) + 1
    cellmax = max(dim1, dim2, dim3)
    # print "dim1 dim2 dim3",dim1,dim2,dim3
    # mat_cvn.make_supercell([dim1,dim2,dim3])
    mat_pos = Poscar(mat_cvn)
    mat_pos.comment = str(p.comment)
    print(mat_pos)
    try:
        toten, final_str, forces = run_job(mat=mat_pos, parameters=parameters)
    except:
        pass
    print("p.comment issssss", p.comment)
    vac = vac_antisite_def_struct_gen(c_size=c_size, struct=final_str)
    def_list, header_list = def_energy(vac=vac, parameters=parameters)
    print(def_list, header_list)
    try:
        print("p.comment issssss", p.comment)
        vac = vac_antisite_def_struct_gen(c_size=c_size, struct=final_str)
        def_list, header_list = def_energy(vac=vac, parameters=parameters)
        print(def_list, header_list)
    except:
        pass
    try:
        surf = pmg_surfer(mat=final_str, min_slab_size=c_size, vacuum=35, max_index=3)
        surf_list, surf_header_list = surf_energy(surf=surf, parameters=parameters)
        print(surf_list, surf_header_list)
    except:
        pass
    cwd = str(os.getcwd())
    if not os.path.exists("Phonon"):
        os.mkdir("Phonon")
    os.chdir("Phonon")
    do_phonons(strt=final_str, parameters=parameters)
    os.chdir(cwd)
    try:
        cwd = str(os.getcwd())
        if not os.path.exists("Phonon"):
            os.mkdir("Phonon")
        os.chdir("Phonon")
        do_phonons(strt=final_str, parameters=parameters)
        os.chdir(cwd)
    except:
        pass
    sub_files = []
    calc = 0
    for a in glob.glob("*.json"):
        fold = a.split(".json")[0]
        cwd = os.getcwd()
        target_file = str(cwd) + str("/") + str(fold)
        dest_file = str(cwd) + str("/") + str(fold) + str(".zip")
        sub_files.append(dest_file)
        ZipDir(
            target_file,
            dest_file,
            contents=[
                "init.mod",
                "potential.mod",
                "in.elastic",
                "data",
                "log.lammps",
                "restart.equil",
                "data0",
            ],
        )
    calc = calc + 1

    target_file = str(cwd)
    dest_file = str(cwd) + str("/") + str("Calc-") + str(calc) + str(".zip")
    ZipDir(target_file, dest_file, contents=sub_files)
    for a in glob.glob("*.json"):
        fold = a.split(".json")[0]
        cwd = os.getcwd()
        dest_file = str(cwd) + str("/") + str(fold) + str(".zip")
        line = str("rm ") + str(dest_file)
        os.system(line)


def main_func(mat=None, parameters={}):
    """
    Call master job function either using mpid or Poscar object
    """
    main(p=mat, parameters=parameters)


def run_job(mat=None, parameters={}, jobname=""):
    """ 
    Generic  function for running LAMMPS job

    Args:
        mat: Poscar object
        parameters: parameters with LAMMPS input information
        jobname: a user-defined jobname
    Returns:
          en: final enery
          final_str: final structure
          forces: forces on final structure
    """
    # def run_job(mat=None,parameters = {'exec':'/cluster/bin/lmp_ctcms-14439-knc6-2','pair_style':'comb3 polar_on','pair_coeff':None,'atom_style': 'charge' ,'control_file':'/home/kamal/inelast.mod'},jobname=''):
    jobname = str(mat.comment)
    # if jobname.startswith('bulk') or jobname.startswith('sbulk'):
    #   parameters['control_file']='/home/kamal/inelast.mod'
    # else:
    #   parameters['control_file']='/home/kamal/inelast_nobox.mod'

    folder = str(os.getcwd()) + str("/") + str(jobname)
    if not os.path.exists(folder):
        os.makedirs(str(jobname))
    os.chdir(str(jobname))
    print("folder name", folder)
    forces = "na"
    if os.path.isfile("./log.lammps"):
        try:

            (
                en,
                press,
                toten,
                c11,
                c22,
                c33,
                c12,
                c13,
                c23,
                c44,
                c55,
                c66,
                c14,
                c16,
                c24,
                c25,
                c26,
                c34,
                c35,
                c36,
                c45,
                c46,
                c56,
            ) = analyz_loge("./log.lammps")
        except:
            pass
    else:
        # pass
        # line=str("#PBS -l pmem=")+str(mem)+'\n'
        # f.write(line)
        write_lammps_data(structure=mat.structure, file="data")
        write_lammps_in(
            structure=mat.structure,
            lammps_in="init.mod",
            lammps_in1="potential.mod",
            lammps_in2="in.main",
            parameters=parameters,
        )
        if "elast" in parameters["control_file"]:
            import jarvis

            base_path = jarvis.__file__.split("/__init__.py")[0]
            displace_file = os.path.join(
                base_path, "lammps", "module_files", "displace.mod"
            )
            shutil.copy2(displace_file, "./")
        os.system(parameters["exec"])
        initial_str = read_data(data="data", ff="potential.mod")

        ####time.sleep(100)
        (
            en,
            press,
            toten,
            c11,
            c22,
            c33,
            c12,
            c13,
            c23,
            c44,
            c55,
            c66,
            c14,
            c16,
            c24,
            c25,
            c26,
            c34,
            c35,
            c36,
            c45,
            c46,
            c56,
        ) = analyz_loge("./log.lammps")

    print("initial,final sr", os.getcwd())
    initial_str = read_data(data="data", ff="potential.mod")
    final_str = read_data(data="data0", ff="potential.mod")
    try:
        forces = read_dump(data="0.dump", ff="potential.mod")
    except:
        pass
    print("initial,final sr2", os.getcwd())

    data_cal = []
    data_cal.append(
        {
            "jobname": jobname,
            "poscar": mat.as_dict(),
            "initial_pos": initial_str.as_dict(),
            "pair_style": str(parameters["pair_style"]),
            "pair_coeff": str(parameters["pair_coeff"]),
            "final_energy": float(toten),
            "en": en,
            "press": press,
            "final_str": final_str.as_dict(),
        }
    )
    json_file = str(jobname) + str(".json")
    os.chdir("../")
    f_json = open(json_file, "w")
    f_json.write(json.dumps(data_cal, indent=4, cls=MontyEncoder))
    f_json.close()
    return en, final_str, forces


if __name__ == "__main__":
    import jarvis

    base_path = jarvis.__file__.split("/__init__.py")[0]
    pos = os.path.join(
        base_path,
        "lammps",
        "examples",
        "Al03.eam.alloy_nist",
        "bulk@mp-134_fold",
        "POSCAR",
    )
    pair_coeff = os.path.join(base_path, "lammps", "examples", "Al03.eam.alloy")

    control_file = os.path.join(base_path, "lammps", "module_files", "inelast.mod")

    mat = Poscar.from_file(pos)

    parameters = {
        "c_size": 3,
        "exec": "/opt/lammps/bin/lammps  <in.main >out",
        "pair_style": "eam/alloy",
        "pair_coeff": pair_coeff,
        "atom_style": "charge",
        "control_file": control_file,
    }
    run_job(mat=mat, parameters=parameters)
    main_func(mat=mat, parameters=parameters)