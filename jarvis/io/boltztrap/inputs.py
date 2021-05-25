"""Class for writing inputs for BoltzTrap calculations."""
import numpy as np
from jarvis.io.vasp.outputs import Vasprun
from jarvis.analysis.structure.spacegroup import Spacegroup3D

Ry_to_ev = 13.6056980659
Angs_to_Bohr = 1.88973


class WriteInputs(object):
    """Write input files for BoltzTrap."""

    def __init__(
        self,
        vasprun_path="",
        dos_type="HISTO",
        tmax=1300,
        tgrid=50,
        doping=[],
        run_type="BOLTZ",
        symprec=1e-2,
        energy_grid=0.005,
        lpfac=10,
        energy_span_around_fermi=1.5,
        energy=None,
        struct=None,
        intrans=None,
    ):
        """
        Require following information.

        energy: energy window.

        struct: Atoms object.

        intrans: name of intrans.

        vasprun_path:  path of vasprun file.

        dos_type: type of densit of states.

        tmax: maximum temperature.

        tgrid: temperature grid.

        doping: doping levels

        run_type:

        symprec: symmetry precision.

        energy_grid: energy grid.

        lpfac:

        energy_span_around_fermi:

        """
        self.energy = energy
        self.struct = struct
        self.intrans = intrans
        self.vasprun_path = vasprun_path
        self.vrun = Vasprun(filename=vasprun_path)
        self.energy_grid = energy_grid
        self.lpfac = lpfac
        self.run_type = run_type
        self.energy_span_around_fermi = energy_span_around_fermi
        self.tmax = tmax
        self.tgrid = tgrid
        self.dos_type = dos_type
        self.doping = doping
        if self.doping == []:
            for d in [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]:
                self.doping.extend([1 * d, 2.5 * d, 5 * d, 7.5 * d])
            self.doping.append(1e22)

    def write_intrans(self, filename="boltztrap.intrans"):
        """Write BoltzTrap input intrans file."""
        f = open(filename, "w")
        scissor = 0.0
        nelect = int(float(self.vrun.all_input_parameters["NELECT"]))
        setgap = 1 if scissor > 0.0001 else 0
        f.write("GENE          # use generic interface\n")
        f.write(
            "1 0 %d %f         # iskip (not presently used) idebug "
            "setgap shiftgap \n" % (setgap, scissor / float(Ry_to_ev))
        )
        f.write(
            "0.0 %f %f %6.1f     # Fermilevel (Ry),energygrid,energy "
            "span around Fermilevel, number of electrons\n"
            % (
                self.energy_grid / float(Ry_to_ev),
                self.energy_span_around_fermi / float(Ry_to_ev),
                nelect,
            )
        )
        f.write(
            "CALC                    # CALC (calculate expansion "
            "coeff), NOCALC read from file\n"
        )
        f.write(
            "%d                        # lpfac, number of latt-points "
            "per k-point\n" % self.lpfac
        )
        f.write(
            "%s                     # run mode (only BOLTZ is "
            "supported)\n" % self.run_type
        )

        f.write(
            ".15                       # (efcut) energy range of "
            "chemical potential\n"
        )
        f.write(
            "{} {}                  # Tmax, temperature grid\n".format(
                self.tmax, self.tgrid
            )
        )

        f.write(
            "-1.  # energyrange of bands given DOS output sig_xxx and "
            "dos_xxx (xxx is band number)\n"
        )

        f.write(self.dos_type + "\n")
        self.taurf = 0
        self.tauexp = 0
        self.tauen = 0
        f.write("{} {} {} 0 0 0\n".format(self.taurf, self.tauexp, self.tauen))

        f.write("{}\n".format(2 * len(self.doping)))
        for i in self.doping:
            f.write(str(i) + "\n")
        for i in self.doping:
            f.write(str(-i) + "\n")
        f.close()

    def write_struct(self, filename="boltztrap.struct"):
        """Write BoltzTrap based struct file."""
        atoms = self.vrun.all_structures[-1]
        spg = Spacegroup3D(atoms)
        spg_symb = spg.space_group_symbol
        formula = atoms.composition.formula
        operations = spg._dataset["rotations"]
        lattice_mat = np.array(atoms.lattice_mat) * Angs_to_Bohr
        f = open(filename, "w")
        f.write("%s %s\n" % (formula, spg_symb))
        f.write(
            "%12.8f %12.8f %12.8f\n"
            % (lattice_mat[0][0], lattice_mat[0][1], lattice_mat[0][2])
        )
        f.write(
            "%12.8f %12.8f %12.8f\n"
            % (lattice_mat[1][0], lattice_mat[1][1], lattice_mat[1][2])
        )
        f.write(
            "%12.8f %12.8f %12.8f\n"
            % (lattice_mat[2][0], lattice_mat[2][1], lattice_mat[2][2])
        )
        f.write("%d\n" % (len(operations)))
        for c in operations:
            for row in c:
                f.write("{}\n".format(" ".join(str(i) for i in row)))
        f.close()

    def write_energy(self, filename="boltztrap.energyso", trim=0.1):
        """Write energy information from DFT."""
        kpoints = self.vrun.kpoints._kpoints
        eigs_up, eigs_dn = self.vrun.eigenvalues
        ef = self.vrun.efermi
        target = 2 * int(len(eigs_dn[0]) * (1 - trim))  # +1
        # print("target", target)
        f = open(filename, "w")
        line = str("system \n") + str(len(kpoints)) + "\n"
        f.write(line)
        for i, j, k in zip(kpoints, eigs_up, eigs_dn):
            count = 0
            line = (" ".join(map(str, i))) + str(" ") + str(target) + "\n"
            # f.write(line)
            f.write("%12.8f %12.8f %12.8f %d\n" % (i[0], i[1], i[2], target))
            for m, n in zip(j, k):
                count = count + 2
                if count <= target:
                    en_up = round((m[0] - ef) / float(Ry_to_ev), 8)
                    en_dn = round((n[0] - ef) / float(Ry_to_ev), 8)
                    f.write("%18.8f\n" % (en_up))
                    f.write("%18.8f\n" % (en_dn))
        f.close()


"""
if __name__ == "__main__":
    from jarvis.io.vasp.outputs import Vasprun

    vrun = Vasprun(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/vasprun.xml"
    )
    print(vrun.final_energy)
    inp = WriteInputs(
        vasprun_path="/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/vasprun.xml"
    )
    inp.write_energy()
    inp.write_struct()
    inp.write_intrans()
    cmd = "~/anaconda2/bin/x_trans BoltzTraP -so"
    os.system(cmd)
"""
