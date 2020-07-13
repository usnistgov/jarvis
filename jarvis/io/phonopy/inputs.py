"""Modules to generate input files for phonopy using Atoms object."""

import numpy as np
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.kpoints import Kpoints3D


class PhonopyInputs(object):
    """
    Generate input files for phonopy post-processing.

    Make sure to obtain FORCE_CONSTANTS first
    E.g. run phonopy --fc vasprun.xml first.
    """

    def __init__(self, atoms):
        """Intialize with jarvis.core.Atoms class."""
        self.atoms = atoms
        self.tag = " ".join(list(set(atoms.elements)))

    def prim_axis(self):
        """
        Use when taking conventional cell.

        Such as elastic constant calculations.
        """
        prim = self.atoms.get_primitive_atoms
        if prim.num_atoms < self.atoms.num_atoms:
            spaceg = Spacegroup3D(self.atoms)
            latt = spaceg.lattice_system
            sgp = spaceg.space_group_symbol
            if latt == "rhombohedral":
                transf = (
                    np.array(
                        [[2, -1, -1], [1, 1, -2], [1, 1, 1]], dtype=np.float
                    )
                    / 3
                )
            elif "I" in sgp:
                transf = (
                    np.array(
                        [[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.float
                    )
                    / 2
                )

            elif "F" in sgp:
                transf = (
                    np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=np.float)
                    / 2
                )
            elif "A" in sgp:
                transf = (
                    np.array(
                        [[2, 0, 0], [0, 1, -1], [0, 1, 1]], dtype=np.float
                    )
                    / 2
                )
            elif "C" in sgp:
                if latt == "monoclinic":
                    transf = (
                        np.array(
                            [[1, 1, 0], [-1, 1, 0], [0, 0, 2]], dtype=np.float
                        )
                        / 2
                    )
                else:
                    transf = (
                        np.array(
                            [[1, -1, 0], [1, 1, 0], [0, 0, 2]], dtype=np.float
                        )
                        / 2
                    )
        else:
            transf = np.eye(3)
        ax = ""
        for el in transf:
            ax = (
                str(ax)
                + str(" ")
                + str(el[0])
                + str(" ")
                + str(el[1])
                + str(" ")
                + str(el[2])
            )
        ax_line = str("PRIMITIVE_AXIS = ") + str(ax) + "\n"
        return ax_line

    def mesh_dos(
        self,
        filename="meshdos.conf",
        dim=[1, 1, 1],
        factor=521.471,
        grid=[31, 31, 31],
    ):
        """
        Use for DOS and thermal properties.

        phonopy -p meshdos.conf
        phonopy -t meshdos.conf
        """
        mesh = open(filename, "w")
        line = str("FORCE_CONSTANTS = READ") + "\n"
        mesh.write(line)
        line = str("FREQUENCY_CONVERSION_FACTOR = ") + str(factor) + "\n"
        mesh.write(line)
        line = str("ATOM_NAME = ") + str(self.tag) + "\n"
        mesh.write(line)
        line = str("DIM = ") + " ".join(map(str, dim)) + "\n"
        mesh.write(line)
        line = str("MP = ") + " ".join(map(str, grid)) + "\n"
        mesh.write(line)
        mesh.close()

    def mesh_bands(
        self,
        dim=[1, 1, 1],
        filename="band.conf",
        factor=521.471,
        line_density=20,
    ):
        """
        Use for making phonon bandstructure plot.

        After running phonopy -p bandd.conf, bandplot -o PBAND.png
        """
        kpoints = Kpoints3D().kpath(self.atoms, line_density=line_density)
        all_kp = kpoints._kpoints
        labels = kpoints._labels
        all_labels = ""
        all_lines = ""
        for lb in labels:
            if lb == "":
                lb = None
            all_labels = all_labels + str(lb) + str(" ")
        for k in all_kp:
            all_lines = (
                all_lines
                + str(k[0])
                + str(" ")
                + str(k[1])
                + str(" ")
                + str(k[2])
                + str(" ")
            )
        file = open(filename, "w")
        line = str("FREQUENCY_CONVERSION_FACTOR = ") + str(factor) + "\n"
        file.write(line)
        ax_line = self.prim_axis()
        file.write(ax_line)
        line = str("ATOM_NAME = ") + str(self.tag) + "\n"
        file.write(line)
        line = str("DIM = ") + " ".join(map(str, dim)) + "\n"
        file.write(line)
        line = str("FORCE_CONSTANTS = READ") + "\n"
        file.write(line)
        line = str("BAND= ") + str(all_lines) + "\n"
        file.write(line)
        line = str("BAND_LABELS= ") + str(all_labels) + "\n"
        file.write(line)
        file.close()

    def mesh_irreps(
        self, factor=521.471, filename="irreps.conf", dim=[1, 1, 1], tol=1e-3
    ):
        """Use for irreducible phonon representation."""
        file = open(filename, "w")
        line = str("IRREPS = 0 0 0 ") + str(tol) + "\n"
        file.write(line)
        line = str("FREQUENCY_CONVERSION_FACTOR = ") + str(factor) + "\n"
        file.write(line)
        ax_line = self.prim_axis()
        file.write(ax_line)
        line = str("SHOW_IRREPS = .TRUE.") + "\n"
        file.write(line)
        line = str("DIM = ") + " ".join(map(str, dim)) + "\n"
        file.write(line)
        line = str("FORCE_CONSTANTS = READ") + "\n"
        file.write(line)
        file.close()

    def animation(self, factor=521.471, filename="anim.conf", dim=[1, 1, 1]):
        """Use for atom dancing/movement animation."""
        anim = open(filename, "w")
        line = str("FREQUENCY_CONVERSION_FACTOR = ") + str(factor) + "\n"
        anim.write(line)
        line = str("FORCE_CONSTANTS = READ") + "\n"
        anim.write(line)
        line = str("ATOM_NAME = ") + str(self.tag) + "\n"
        anim.write(line)
        line = str("ANIME_TYPE = JMOL") + "\n"
        anim.write(line)
        line = str("ANIME = 0 7 0  0 0 0") + "\n"
        anim.write(line)
        line = str("DIM = ") + " ".join(map(str, dim)) + "\n"
        anim.write(line)
        anim.close()

    def generate_all_files(self):
        """Generate all the input files mentioned above."""
        self.mesh_dos()
        self.mesh_bands()
        self.mesh_irreps()
        self.animation()


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    PhonopyInputs(atoms=Si).generate_all_files()
"""
