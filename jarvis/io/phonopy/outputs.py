"""Module for post-processing phonopy outputs."""
import yaml
import matplotlib.pyplot as plt
from yaml import Loader
import numpy as np
from jarvis.io.phonopy.inputs import PhonopyInputs

# from jarvis.analysis.structure.spacegroup import Spacegroup3D
# from jarvis.io.wannier.outputs import WannierHam
from jarvis.io.phonopy.fcmat2hr import get_phonon_hr

# VaspToTHz = 15.633302300230191
try:
    from phonopy import Phonopy
    from phonopy import load
    from phonopy.structure.cells import determinant
    from phonopy.structure.cells import get_reduced_bases
except Exception as exp:
    print("Phonopy is not installed.", exp)
    pass


from math import pi as pi

"""
Constants
"""
kB = 1.38064852e-23
hbar = 1.0545718e-34
h = 6.62607004e-34
Na = 6.0221409e23
e = 1.60218e-19


def bandstructure_plot(band_yaml="", plot=False):
    """Get phonopy bandstructure info."""
    with open(band_yaml, "r") as f:
        data = yaml.load(f, Loader=Loader)
    frequencies = []
    distances = []
    qpoints = []
    label_points = []
    labels = []
    for j, v in enumerate(data["phonon"]):
        if "label" in v and v["label"] != "None":
            labels.append(v["label"])
            label_points.append(v["distance"])
        frequencies.append([f["frequency"] for f in v["band"]])
        qpoints.append(v["q-position"])
        distances.append(v["distance"])
    if plot:
        for i in range(np.array(frequencies).shape[1]):
            plt.plot(distances, np.array(frequencies)[:, i])
        plt.xticks(label_points, labels)
    return frequencies, distances, labels, label_points


def total_dos(tot_dos="", plot=False):
    """Get total dos info."""
    f = open(tot_dos, "r")
    freq = []
    pdos = []
    for lines in f.readlines():
        if not str(lines.split()[0]).startswith("#"):
            #   print (lines)
            # else:
            freq.append(float(lines.split()[0]))
            pdos.append(float(lines.split()[1]))
    if plot:
        plt.plot(freq, pdos)
    return freq, pdos


"""
More generalized read_fcmethod. But need to fix.
"""


def read_fc(filename="FORCE_CONSTANTS"):
    """Read phonopy generated force constants."""
    f = open(filename, "r")
    lines = f.read().splitlines()
    f.close()
    n_patoms = int(lines[0].split()[0])
    try:
        n_satoms = int(lines[0].split()[1])
    except BaseException:
        n_satoms = n_patoms
    else:
        print("Error fetching supercell atoms.")
    fc = np.zeros((n_patoms, n_satoms, 3, 3), dtype="double")
    # print ('natoms=',natoms)
    patom_id = 0
    satom_id = 0
    for ii, i in enumerate(lines):
        if ii > 0 and ii % 4 == 0:
            atoms_ids = [int(a) for a in lines[ii - 3].split()]
            print(atoms_ids)
            vals = str(lines[ii - 2]) + " " + str(lines[ii - 1])
            vals = vals + " " + (lines[ii])
            vals = np.array(vals.split(), dtype="double").reshape(3, 3)
            fc[patom_id, satom_id] = vals
            satom_id += 1
            if satom_id == n_satoms:
                satom_id = 0
                patom_id += 1
    #            fc[atoms_ids[0] - 1, atoms_ids[1] - 1] = vals
    return fc


def get_Phonopy_obj(
    atoms,
    phonopy_yaml=None,
    FC_file=None,
    unitcell_file=None,
    scell_file=None,
    factor=None,
    symprec=1e-05,
    scell=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
):
    """
    Create separate method for generating a Phonopy object
    """
    if phonopy_yaml is not None:
        phonon = load(phonopy_yaml, force_constants_filename=FC_file)
    elif scell_file is not None:
        unitcell = atoms.phonopy_converter()
        prim_mat = np.array(
            PhonopyInputs(atoms).prim_axis().split("=")[1].split(), 
		dtype="float"
        ).reshape(3, 3)
        if factor is None:
            from phonopy.units import VaspToCm

            factor = VaspToCm
        phonon = load(
            unitcell=unitcell,
            # unitcell_filename = unitcell_file,
            # supercell_matrix = scell,
            supercell_filename=scell_file,
            force_constants_filename=FC_file,
            symprec=symprec,
            is_symmetry=True,
            log_level=2,
        )
    else:
        raise Exception("Invalid inputs to generate Phonopy object.")
    try:
        dim = np.array([scell[i][i] for i in range(3)])
        setattr(phonon, "scell_dim", dim)
    except BaseException:
        setattr(phonon, "scell_dim", None)
    else:
        print("Error setting supercell dimensions.")
    return phonon


def get_thermal_properties(phonon_obj,
	mesh=[1, 1, 1],
	tmin=0,
	tmax=100,
	step=10):
    """
    Returns dictionary of thermal properties, including Helmholtz free energy,
    vibrational entropy, and heat capacity

    Parameters
    ----------
    phonon_obj : Phonopy
        phonopy object
    mesh : list, optional
        Mesh size for thermal property calculation. The default is [1, 1, 1].
    tmin : float, optional
        Minimum temperature. The default is 0.
    tmax : float, optional
        Maximum temperature. The default is 100.
    step : float, optional
        Temperature step size. The default is 10.

    Returns
    -------
    tp_dict : TYPE
        DESCRIPTION.

    """
    phonon_obj.run_mesh(mesh)
    phonon_obj.run_thermal_properties(t_step=step, t_max=tmax, t_min=tmin)
    tp_dict = phonon_obj.get_thermal_properties_dict()
    return tp_dict


def get_modal_heat_capacity(phonon_obj, mesh=[1, 1, 1], T=300):
    phonon_obj.run_mesh(mesh)
    mesh_dict = phonon_obj.get_mesh_dict()
    omega = np.array(mesh_dict["frequencies"]) * 1e12
    x = h * omega / (kB * T)  # omega is ordinal not angular
    mode_C = (kB / e) * (x) ** 2 * (np.exp(x) / ((np.exp(x) - 1) ** 2))
    return mode_C


def get_phonon_tb(
    # phonopy_atoms=[],
    atoms=[],
    fc=[],
    out_file="phonopyTB_hr.dat",
    distance_to_A=1.0,
    scell=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    factor=None,
    symprec=1e-05,
    displacement_distance=0.01,
):
    """Generate phonon TB Hamiltonia, along with WannierHamn."""
    if factor is None:
        from phonopy.units import VaspToCm

        factor = VaspToCm
    # Forked from Wannier-tools
    unitcell = atoms.phonopy_converter()
    # unitcell = phonopy_atoms
    prim_mat = np.array(
        PhonopyInputs(atoms).prim_axis().split("=")[1].split(), dtype="float"
    ).reshape(3, 3)
    # print("cell", unitcell.cell)
    num_atom = unitcell.get_number_of_atoms()
    num_satom = determinant(scell) * num_atom
    if fc.shape[0] != num_satom:
        print("Check Force constant matrix.")

    phonon = Phonopy(
        unitcell,
        scell,
        primitive_matrix=prim_mat,
        factor=factor,
        dynamical_matrix_decimals=None,
        force_constants_decimals=None,
        symprec=symprec,
        is_symmetry=True,
        use_lapack_solver=False,
        log_level=1,
    )
    print("Phonopy object generated")
    supercell = phonon.get_supercell()
    primitive = phonon.get_primitive()
    # Set force constants
    phonon.set_force_constants(fc)
    phonon._set_dynamical_matrix()
    print("here")
    dmat = phonon._dynamical_matrix
    print(dmat)
    # rescale fcmat by THZ**2
    fcmat = dmat._force_constants * factor ** 2  # FORCE_CONSTANTS
    # fcmat = dmat._force_constants * factor ** 2  # FORCE_CONSTANTS
    smallest_vectors = dmat._smallest_vectors
    # mass = dmat._mass
    mass = dmat._pcell.get_masses()
    print("mass=", mass)
    multi = dmat._multiplicity
    reduced_bases = get_reduced_bases(supercell.get_cell(), symprec)
    positions = np.dot(supercell.get_positions(), np.linalg.inv(reduced_bases))
    # for pos in positions: pos -= np.rint(pos)
    relative_scale = np.dot(reduced_bases, np.linalg.inv(primitive.get_cell()))
    super_pos = np.zeros((num_satom, 3), dtype=np.float64)
    for i in range(num_satom):
        super_pos[i] = np.dot(positions[i], relative_scale)
    p2s_map = dmat._p2s_map = primitive.get_primitive_to_supercell_map()
    s2p_map = dmat._s2p_map = primitive.get_supercell_to_primitive_map()
    num_satom = supercell.get_number_of_atoms()
    num_patom = primitive.get_number_of_atoms()
    get_phonon_hr(
        fcmat,
        smallest_vectors,
        mass,
        multi,
        super_pos,
        p2s_map,
        s2p_map,
        num_satom,
        num_patom,
        out_file,
    )
    print("phonopy_TB.dat generated! ")


if __name__ == "__main__":
    from phonopy.interface.vasp import read_vasp
    from jarvis.core.atoms import Atoms

    test_dir = "Si-testing/"
    pos = test_dir + "POSCAR-unitcell"
    fc_file = test_dir + "FORCE_CONSTANTS"
    a = Atoms.from_poscar(pos)
    fc = read_fc(fc_file)
    phonopy_atoms = read_vasp(pos)
    phonon_obj = get_Phonopy_obj(
        a,
        phonopy_yaml="Si-testing/phonopy.yaml",
        FC_file=fc_file,
        scell=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    )
    phonon_obj.run_mesh(
        [11, 11, 11], with_eigenvectors=True, with_group_velocities=True
    )
    mesh_dict = phonon_obj.get_mesh_dict()
    phonon_obj.run_total_dos()
    phonon_obj.plot_total_dos().show()
    tp_dict = get_thermal_properties(
        phonon_obj, mesh=[11, 11, 11], tmin=0, tmax=300, step=100
    )
#    gv_by_gv = get_gv_outer_product(phonon_obj, mesh=[11, 11, 11])
#    get_phonon_tb(fc=fc, atoms=a)
#    cvn = Spacegroup3D(a).conventional_standard_structure
#    w = WannierHam("phonopyTB_hr.dat")
#    w.get_bandstructure_plot(atoms=cvn, yrange=[0, 550])
