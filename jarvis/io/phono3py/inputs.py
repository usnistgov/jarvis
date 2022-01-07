import numpy as np
import subprocess
import spglib
from jarvis.io.phonopy.outputs import get_Phonopy_obj
import os

#Maybe add get_gridpoints function here?

#Should be a way to get supercell dimensions programmatically?
def prepare_jdos(phonopy_obj, poscar = "", mesh = [1,1,1],\
                 scell_dim = [1,1,1], temperatures = None,\
                     num_freq_points = 201, pa = None, run = False):
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
    mapping, grid = spglib.get_ir_reciprocal_mesh(
        mesh, cell, is_shift=[0, 0, 0]
    )
    gridpt_uids = np.unique(mapping)
    dim = phonopy_obj.scell_dim
    phono3py_cmd =     phono3py_cmd = 'phono3py --fc2 --dim="' +\
        " ".join([str(d) for d in dim]) + '" --mesh="' +\
            " ".join([str(d) for d in mesh]) + '" -c ' + poscar + ' --jdos' +\
    ' --gp="' + ' '.join([str(g) for g in gridpt_uids]) +\
        '" --num-freq-points=' + str(num_freq_points)
    
    if temperatures is not None:
        phono3py_cmd = phono3py_cmd + '" --ts="' + " ".join(
        [str(d) for d in temperatures])
    
    if pa is not None:
        phono3py_cmd = phono3py_cmd + ' --pa= "' +\
            " ".join([str(d) for d in pa]) + '"'
    
    print(phono3py_cmd)
    
    if run:
        subprocess.call(phono3py_cmd, shell=True)
    
    
if __name__ == '__main__':
    from jarvis.core.atoms import Atoms

    test_dir = "Si-testing/"
    os.chdir(test_dir)
    pos = "POSCAR-unitcell"
    atoms = Atoms.from_poscar(pos)
    phonon_obj = get_Phonopy_obj(
        atoms,
        phonopy_yaml="phonopy.yaml",
        FC_file="FORCE_CONSTANTS",
        scell=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    )
    prepare_jdos(phonon_obj, poscar = pos, mesh = [11, 11, 11],\
                 scell_dim = [2,2,2], run = True)
