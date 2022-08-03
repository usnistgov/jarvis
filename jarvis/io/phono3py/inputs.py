import numpy as np
import subprocess
import spglib
from jarvis.io.phonopy.outputs import get_Phonopy_obj
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.io.vasp.inputs import Poscar
import os
from phono3py import Phono3pyJointDos

# Maybe add get_gridpoints function here?

# Should be a way to get supercell dimensions programmatically?
def prepare_jdos(
    phonopy_obj,
    poscar="",
    mesh=[1, 1, 1],
    scell_dim=None,
    temperatures=None,
    num_freq_points=201,
    pa=None,
    run=False,
):
    """
    

    Parameters
    ----------
    phonopy_obj : TYPE
        DESCRIPTION.
    poscar : TYPE, optional
        DESCRIPTION. The default is "".
    scell_dim : TYPE, optional
        DESCRIPTION. The default is [1,1,1].
    temperatures : TYPE, optional
        DESCRIPTION. The default is None.
    num_freq_points : TYPE, optional
        DESCRIPTION. The default is 201.
    pa : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    latt_vecs = phonopy_obj.get_primitive().get_cell()
    positions = phonopy_obj.get_primitive().get_scaled_positions()
    atom_type = phonopy_obj.get_primitive().get_atomic_numbers()
    cell = (latt_vecs, positions, atom_type)
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
    gridpt_uids = np.unique(mapping)
    dim = phonopy_obj.scell_dim
    phono3py_cmd = (
        'phono3py --fc2 --dim="'
        + " ".join([str(d) for d in dim])
        + '" --mesh="'
        + " ".join([str(d) for d in mesh])
        + '" -c '
        + poscar
        + " --jdos"
        + ' --gp="'
        + " ".join([str(g) for g in gridpt_uids])
        + '" --num-freq-points='
        + str(num_freq_points)
    )

    if temperatures is not None:
        phono3py_cmd = (
            phono3py_cmd + '" --ts="' + " ".join([str(d) for d in temperatures])
        )

    if pa is not None:
        phono3py_cmd = phono3py_cmd + ' --pa= "' + " ".join([str(d) for d in pa]) + '"'

    print(phono3py_cmd)

    if run:
        subprocess.call(phono3py_cmd, shell=True)


# Try to use the Phono3pyJDOS class?
def prepare_jdos_api(phonopy_obj, mesh = [1, 1, 1], temperatures=None,\
    num_freq_points=201, pa=None, nac_params = None,
    write_jdos = True):
    '''
    Run JDOS calculation using the Phono3pyJointDOS Class
    Currently assumes phonopy_obj contains a force cosntants attribute

    Parameters
    ----------
    phonopy_obj : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    ph3_jdos = Phono3pyJointDos(phonopy_obj._supercell,\
                                phonopy_obj._primitive,\
                                    phonopy_obj._force_constants,\
                                        mesh = mesh,\
                                            temperatures = temperatures,\
                                                nac_params = nac_params,\
                                                    num_frequency_points = num_freq_points,\
                                                        log_level = 2)
#    latt_vecs = phonopy_obj.get_primitive().get_cell()
#    positions = phonopy_obj.get_primitive().get_scaled_positions()
#    atom_type = phonopy_obj.get_primitive().get_atomic_numbers()
#    cell = (latt_vecs, positions, atom_type)
#    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
#    gridpt_uids = np.unique(mapping)
    ga, gp, mapping = phonopy_obj.get_mesh_grid_info()    
    ph3_jdos.run(gp, write_jdos = write_jdos)
    


# Should add a prepare_gruneisen_quasiharmonic function as well


def prepare_gruneisen_quasiharmonic(orig_poscar: str, scale):
    """
    Generates the dilated and constricted POSCAR files for a VASP run
    scale : should be larger than 1, scales to the dilated unit cell
    """
    # Read original POSCAR
    orig_file = open(orig_poscar, "r")
    pos_text = orig_file.read().splitlines()
    scale_plus = float(pos_text[1]) * scale
    scale_minus = float(pos_text[1]) / scale
    # First write the large unit cell POSCAR
    pos_text[1] = str(scale_plus)
    new_file = open("POSCAR-plus", "w")
    new_file.writelines("\n".join(pos_text))
    # Next write small unit cell POSCAR
    pos_text[1] = str(scale_minus)
    new_file = open("POSCAR-minus", "w")
    new_file.writelines("\n".join(pos_text))


def prepare_gruneisen_FC3(
    phonopy_obj,
    poscar="",
    band_calc=False,
    line_density=1,
    mesh=[1, 1, 1],
    nac=False,
    run=False,
    plot=False,
):
    """
    

    Parameters
    ----------
    phonopy_obj : TYPE
        DESCRIPTION.
    poscar : TYPE, optional
        DESCRIPTION. The default is "".
    band_calc : boolean, optional
    If true, performs a band calculation instead of mesh calculation. Default
    is false.
    mesh : TYPE, optional
        DESCRIPTION. The default is [1,1,1].

    Returns
    -------
    None.

    """
    dim = phonopy_obj.scell_dim
    phono3py_cmd = (
        'phono3py --fc3 --fc2 --dim="'
        + " ".join([str(d) for d in dim])
        + '" -v -c '
        + poscar
        + " --gruneisen"
    )

    if nac:
        phono3py_cmd = phono3py_cmd + " --nac"

    if band_calc:
        atoms = Atoms.from_poscar(poscar)
        kpoints = Kpoints3D().kpath(atoms, line_density=line_density)
        all_kp = kpoints._kpoints
        all_lines = ""
        for k in all_kp:
            all_lines = (
                all_lines
                + str(k[0])
                + str(" ")
                + str(k[1])
                + str(" ")
                + str(k[2])
                + str("  ")
            )
        phono3py_cmd = phono3py_cmd + ' --band="' + all_lines + '"'
    else:
        phono3py_cmd = (
            phono3py_cmd + ' --mesh="' + " ".join([str(d) for d in mesh])
        ) + '"'
    print(phono3py_cmd)
    if run:
        subprocess.call(phono3py_cmd, shell=True)


if __name__ == "__main__":

    test_dir = "Si-testing/phono3py-example-Si-PBEsol/"
    os.chdir(test_dir)
    pos = "POSCAR-unitcell"
    atoms = Atoms.from_poscar(pos)
    phonon_obj = get_Phonopy_obj(
        atoms,
        phonopy_yaml="phono3py.yaml",
        FC_file="fc2.hdf5",
        scell=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    )
    phonon_obj.run_mesh([11, 11, 11], with_group_velocities=True)
    mesh_dict = phonon_obj.get_mesh_dict()
    # prepare_jdos(
    #     phonon_obj, poscar=pos, mesh=[11, 11, 11], scell_dim=[2, 2, 2], run=True
    # )
    prepare_jdos_api(phonon_obj, mesh=[11, 11, 11], write_jdos = True)
    # prepare_gruneisen_FC3(
    #     phonon_obj, poscar=pos, mesh=[2, 2, 2], band_calc=False, run=True, plot=True
    # )
    # prepare_gruneisen_quasiharmonic("POSCAR-unitcell", 1.00335)
