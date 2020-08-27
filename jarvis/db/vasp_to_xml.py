"""Module to make XML for JARVIS-REST-API schema."""

import os, glob
from jarvis.core.atoms import Atoms
import yaml
from jarvis.core.atoms import get_supercell_dims
from jarvis.db.jsonutils import loadjson
import io
from jarvis.io.vasp.outputs import Oszicar, Vasprun, Outcar
from matplotlib.pyplot import imread
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.analysis.topological.spillage import Spillage
from jarvis.analysis.phonon.ir import ir_intensity
from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.analysis.thermodynamics.energetics import form_enp
from jarvis.core.utils import recast_array_on_uniq_array_elements
import numpy as np
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.io.boltztrap.outputs import BoltzTrapOutput
from jarvis.core.utils import stringdict_to_xml, array_to_string


class VaspToApiXmlSchema(object):
    def __init__(self, folder="", data={}):
        self.folder = folder
        self.data = data

    def ebandstruct(self, vrun="", kp=""):
        line = ""
        try:
            f = open(kp, "r")
            lines = f.read().splitlines()
            f.close()

            vrun = Vasprun(vrun)

            kp_labels = []
            kp_labels_points = []
            for ii, i in enumerate(lines):
                if ii > 2:
                    tmp = i.split()
                    if len(tmp) == 5:
                        tmp = str("$") + str(tmp[4]) + str("$")
                        if len(kp_labels) == 0:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)
                        elif tmp != kp_labels[-1]:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)
            tmp = 0.0
            zero_efermi = True
            if zero_efermi:
                tmp = float(vrun.efermi)
            spin = 0
            up_bands_y = []
            up_bands_x = []
            for i, ii in enumerate(vrun.eigenvalues[spin][:, :, 0].T - tmp):
                y = ",".join(map(str, ii))
                x = ",".join(map(str, range(0, len(ii))))
                up_bands_y.append(y)
                up_bands_x.append(x)
            spin = 1
            down_bands_y = []
            down_bands_x = []
            for i, ii in enumerate(vrun.eigenvalues[spin][:, :, 0].T - tmp):
                y = ",".join(map(str, ii))
                x = ",".join(map(str, range(0, len(ii))))
                down_bands_y.append(y)
                down_bands_x.append(x)
            line = ""
            line += (
                '<kp_labels_points>"'
                + ",".join(map(str, kp_labels_points))
                + '"</kp_labels_points>'
            )
            line += (
                '<kp_labels>"'
                + ",".join(map(str, kp_labels))
                + '"</kp_labels>'
            )
            line += (
                '<spin_up_bands_x>"'
                + ";".join(map(str, up_bands_x))
                + '"</spin_up_bands_x>'
            )
            line += (
                '<spin_up_bands_y>"'
                + ";".join(map(str, up_bands_y))
                + '"</spin_up_bands_y>'
            )
            line += (
                '<spin_down_bands_x>"'
                + ";".join(map(str, down_bands_x))
                + '"</spin_down_bands_x>'
            )
            line += (
                '<spin_down_bands_y>"'
                + ";".join(map(str, down_bands_y))
                + '"</spin_down_bands_y>'
            )
        except Exception:
            print("Could not get bandstructure info,vrun,kp", vrun, kp)
            pass

        return line

    def encut_kp(self):
        folder = self.folder
        os.chdir(folder)
        info = ""

        try:
            encut_files = []
            encut_values = []
            encut_based_energies = []
            for i in glob.glob("*.json"):
                if "ENCUT" in i:
                    encut_values.append(
                        int(str(i.split("-")[-1]).split(".json")[0])
                    )
                    encut_files.append(i)
                    erun = os.path.join(
                        folder, i.split(".json")[0], "vasprun.xml"
                    )
                    energy = Vasprun(erun).final_energy
                    encut_based_energies.append(energy)

            kplength_files = []
            kp_values = []
            kp_based_energies = []
            for i in glob.glob("*.json"):
                if "KPOINT" in i:
                    kp_values.append(
                        int(str(i.split("-")[-1]).split(".json")[0])
                    )
                    kplength_files.append(i)
                    krun = os.path.join(
                        folder, i.split(".json")[0], "vasprun.xml"
                    )
                    energy = Vasprun(krun).final_energy
                    kp_based_energies.append(energy)

            info = {}
            info["kp_values"] = ",".join(map(str, kp_values))
            info["kp_based_energies"] = ",".join(map(str, kp_based_energies))
            info["encut_values"] = ",".join(map(str, encut_values))
            info["encut_based_energies"] = ",".join(
                map(str, encut_based_energies)
            )
            self.encut_kp_info = info
        except Exception:
            print("Check the convergence folders are named properly")
            pass
        os.chdir(folder)

        return info

    def electronic_dos_info(self, vrun):
        """
            Use , ; _ @ # $ & special-character symbols to condens data to strings i.e. 7 leve deep.

            Note in the JARVIS-API XML schema we use strings only to store information.
            """
        # , for individual array element
        # ; same plot
        # @ different plot
        # _ title vs rest
        line = ""
        try:
            total_dos = vrun.total_dos
            energies = total_dos[0]
            spdf_dos = vrun.get_spdf_dos()
            atom_dos = atom_resolved_dos = vrun.get_atom_resolved_dos()
           
            # designed as dosname_dosvaluearray_color_label seperated by semicoln for for each plot. List of plots are seperated by @
            line = ""
            for i, ii in enumerate(total_dos):
                if i == 0:
                    line = (
                        line
                        + "<edos_energies>'"
                        + ",".join(map(str, ii))
                        + "'</edos_energies>"
                    )
                elif i == 1:
                    line += (
                        "<total_edos_up>'"
                        + ",".join(map(str, ii))
                        + "'</total_edos_up>"
                    )
                elif i == 2:
                    line += (
                        "<total_edos_down>'"
                        + ",".join(map(str, ii))
                        + "'</total_edos_down>"
                    )

            line += "<spdf_dos>'"
            for i, j in spdf_dos.items():
                line += (
                    "<"
                    + str(i)
                    + ">'"
                    + ",".join(map(str, j))
                    + "'</"
                    + str(i)
                    + ">"
                )

            line += "'</spdf_dos>"
            line += "<elemental_dos>"
            spin_up_info=atom_dos['spin_up_info']
            spin_down_info=atom_dos['spin_down_info']
            line+='<spin_up_info>"'
            for i,j in spin_up_info.items():
                line+=str(i)+str('_')+",".join(map(str, j))+';'
            line+='"</spin_up_info>'
        
            line+='<spin_down_info>"'
            for i,j in spin_down_info.items():
                line+=str(i)+str('_')+",".join(map(str, j))+';'
            line+='"</spin_down_info>'
            line += "</elemental_dos>"
        except:
            print("Cannot get DOS data", vrun)
        return line

    def electronic_band_struct(self, vasprun="", kpoints_file_path=""):
        line = self.ebandstruct(vrun=vasprun, kp=kpoints_file_path)
        return line

    def main_band(self):
        folder = self.folder
        os.chdir(folder)
        data = ""
        info = {}
        for i in glob.glob("MAIN-BAND*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                kp = os.getcwd() + "/" + folder + "/KPOINTS"
                data = self.electronic_band_struct(
                    vasprun=vrun, kpoints_file_path=kp
                )
            except Exception:
                print("Cannot get ebands info", folder)
                pass
        os.chdir(folder)
        info["main_bands_info"] = data
        return info

    def loptics_optoelectronics(self, vrun=""):
        line = ""
        try:
            reals, imags = Vasprun(vrun).dielectric_loptics
            energies = reals[:, 0]
            line = ""
            line += (
                "<energies>'" + ",".join(map(str, energies)) + "'</energies>"
            )
            for i in np.arange(1, reals.shape[1]):
                line += (
                    "<real_"
                    + str(i)
                    + ">'"
                    + ",".join(map(str, reals[:, i]))
                    + "'</real_"
                    + str(i)
                    + ">"
                )
                line += (
                    "<imag_"
                    + str(i)
                    + ">'"
                    + ",".join(map(str, imags[:, i]))
                    + "'</imag_"
                    + str(i)
                    + ">"
                )
        except Exception:
            print("Cannot get loptics data", vrun)
            pass

        return line

    def main_optics_semilocal(self):
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-OPTICS*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                data = self.loptics_optoelectronics(vrun)
            except Exception:
                print("Cannot get semilocal optics data", folder)
                pass
        os.chdir(folder)
        info["main_optics_info"] = data
        self.main_optics_semilocal_info = data
        return info

    def main_optics_mbj(self):
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-MBJ*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                data = self.loptics_optoelectronics(vrun)
            except:
                print("Cannot get mbj optics data", folder)
                pass
        os.chdir(folder)
        info["main_optics_mbj_info"] = data
        self.main_optics_mbj_info = data
        return info

    def spillage(
        self,
        soc_wf="MAIN-SOCSCFBAND-bulk@JVASP-1067_mp-541837/WAVECAR",
        nonsoc_wf="MAIN-MAGSCFBAND-bulk@JVASP-1067_mp-541837/WAVECAR",
    ):

        line = ""
        try:
            sp = Spillage(wf_noso=nonsoc_wf, wf_so=soc_wf).overlap_so_spinpol()
            soc_bands = self.electronic_band_struct(
                vasprun=soc_wf.replace("WAVECAR", "vasprun.xml"),
                kpoints_file_path=(soc_wf.replace("WAVECAR", "KPOINTS")),
            )
            line += "<soc_bands>" + soc_bands + "</soc_bands>"
            nonsoc_bands = self.electronic_band_struct(
                vasprun=nonsoc_wf.replace("WAVECAR", "vasprun.xml"),
                kpoints_file_path=(nonsoc_wf.replace("WAVECAR", "KPOINTS")),
            )
            line += "<nonsoc_bands>" + nonsoc_bands + "</nonsoc_bands>"
            spillage_k = ",".join(map(str, sp["spillage_k"]))
            line += "<spillage_k>" + spillage_k + "</spillage_k>"
            spillage_kpoints = ",".join(
                [";".join(map(str, i)) for i in sp["kpoints"]]
            )
            line += (
                "<spillage_kpoints>" + spillage_kpoints + "</spillage_kpoints>"
            )
        except Exception:
            print("Cannot get spillage info", soc_wf, nonsoc_wf)
            pass
        return line

    def main_soc_spillage(self):
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-SOCSCFBAND*.json"):

            try:
                folder = i.split(".json")[0]
                soc_wf = os.getcwd() + "/" + folder + "/WAVECAR"
                nonsoc_wf = (
                    os.getcwd()
                    + "/"
                    + folder.replace("SOC", "MAG")
                    + "/WAVECAR"
                )
                data = self.spillage(soc_wf=soc_wf, nonsoc_wf=nonsoc_wf)
            except:
                print("Cannot get spillage info", folder)
                pass
        os.chdir(folder)
        info["main_spillage_info"] = data
        self.main_soc_spillage_info = info
        return info

    def basic_info(
        self,
        id_file="JARVIS-ID",
        material_type="bulk",
        ehull="",
        include_dos_info=True,
        include_neighbor_info=True,
        data_source="JARVIS-DFT-VASP",
    ):
        main_folder = self.folder
        os.chdir(main_folder)
        info = {}
        info["data_source"] = data_source
        info["ehul"] = ehull
        info["source_path"] = main_folder
        info["material_type"] = material_type
        for i in glob.glob("MAIN-RELAX*.json"):
            folder = i.split(".json")[0]
            main_vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
            try:
                vrun = Vasprun(main_vrun)
                atoms = vrun.all_structures[-1]
                f = open(os.path.join(folder, "..", id_file), "r")
                lines = f.read().splitlines()
                f.close()
                id = lines[0]
                info["id"] = id
                formula = atoms.composition.reduced_formula
                info["formula"] = formula
                elements = ",".join(atoms.uniq_species)
                info["elements"] = elements
                method = ""
                try:
                    f = open(os.path.join(folder, "..", "FUNCTIONAL"), "r")
                    lines = f.read().splitlines()
                    f.close()
                    method = lines[0]
                except Exception:
                    print('Cannot find "FUNCTIONAL"')
                    pass
                if method == "PBEBO":
                    method = "OptB88vdW"
                if method == "PBEOR":
                    method = "OptPBEvdW"
                if method == "PBEMK":
                    method = "OptB86bvdW"
                info["method"] = method
                final_energy = vrun.final_energy
                num_atoms = atoms.num_atoms
                fen = ""
                if method == "OptB88vdW":
                    fen = form_enp(atoms=atoms, total_energy=final_energy)
                info["formation_energy"] = fen
                rel_en = final_energy / num_atoms
                info["relaxed_energy"] = round(rel_en, 3)
                spg = Spacegroup3D(atoms)
                prim_atoms = spg.primitive_atoms
                conv_atoms = spg.conventional_standard_structure
                spg_numb = spg.space_group_number
                info["spacegroup_number"] = spg_numb
                info["point_group_symbol"] = spg.point_group_symbol
                spg_symb = spg.space_group_symbol
                info["spg.space_group_symbol"] = spg.space_group_symbol
                crys_system = spg.crystal_system
                info["crys_system"] = crys_system
                conv_params = conv_atoms.lattice.parameters
                a_conv = round(conv_params[0], 2)
                b_conv = round(conv_params[1], 2)
                c_conv = round(conv_params[2], 2)
                info["a_conv"] = a_conv
                info["b_conv"] = b_conv
                info["c_conv"] = c_conv
                alpha_conv = round(conv_params[3], 2)
                beta_conv = round(conv_params[4], 2)
                gamma_conv = round(conv_params[5], 2)
                info["alpha_conv"] = alpha_conv
                info["beta_conv"] = beta_conv
                info["gamma_conv"] = gamma_conv
                prim_params = prim_atoms.lattice.parameters
                a_prim = round(prim_params[0], 2)
                b_prim = round(prim_params[1], 2)
                c_prim = round(prim_params[2], 2)
                info["a_prim"] = a_prim
                info["b_prim"] = b_prim
                info["c_prim"] = c_prim
                alpha_prim = round(prim_params[3], 2)
                beta_prim = round(prim_params[4], 2)
                gamma_prim = round(prim_params[5], 2)
                info["alpha_prim"] = alpha_prim
                info["beta_prim"] = beta_prim
                info["gamma_prim"] = gamma_prim
                prim_natoms = prim_atoms.num_atoms
                conv_natoms = conv_atoms.num_atoms
                info["prim_natoms"] = prim_natoms
                info["conv_natoms"] = conv_natoms
                density = round(atoms.density, 3)
                volume = round(atoms.volume, 3)
                info["density"] = density

                info["volume"] = volume
                packing_fr = ""
                try:
                    packing_fr = atoms.packing_fraction
                except Exception:
                    print("Cannot calculate packing_fr, check atomic_radii")
                    pass
                info["packing_fr"] = packing_fr

                dim = get_supercell_dims(conv_atoms)

                xyz = conv_atoms.make_supercell_matrix(dim).get_xyz_string
                info["xyz"] = '"' + str(xyz).replace("\n", "\\n") + '"'
                if include_neighbor_info:
                    nbr = NeighborsAnalysis(atoms)
                    rdf_bins, rdf_hist, nn = nbr.get_rdf()
                    nbr = NeighborsAnalysis(atoms, max_cut=5.0)
                    ang1_hist, ang1_bins = nbr.ang_dist_first()
                    ang2_hist, ang2_bins = nbr.ang_dist_second()
                    dhd_hist, dhd_bins = nbr.get_ddf()
                    info["rdf_bins"] = "'" + ",".join(map(str, rdf_bins)) + "'"
                    info["rdf_hist"] = "'" + ",".join(map(str, rdf_hist)) + "'"

                    info["ang1_bins"] = ",".join(map(str, ang1_bins))
                    info["ang1_hist"] = ",".join(map(str, ang1_hist))
                    info["ang2_bins"] = ",".join(map(str, ang2_bins))
                    info["ang2_hist"] = ",".join(map(str, ang2_hist))

                scf_indir_gap = vrun.get_indir_gap[0]
                scf_dir_gap = vrun.get_dir_gap
                info["scf_indir_gap"] = round(scf_indir_gap, 3)
                info["scf_dir_gap"] = round(scf_dir_gap, 3)
                try:
                    oszicar = Oszicar(os.path.join(folder, "OSZICAR"))
                    magmom = oszicar.magnetic_moment
                    info["magmom"] = round(float(magmom), 3)
                except:
                    print("Check oszicar")
                    pass
                if include_dos_info:
                    main_relax_dos = self.electronic_dos_info(vrun)
                    info["main_relax_dos"] = main_relax_dos
            except:
                pass
        os.chdir(main_folder)
        self.basic_info_dict = info
        return info

    def dfpt_related(self, vrun="", out=""):
        vrun = Vasprun(vrun)
        data = vrun.dfpt_data
        out = Outcar(out)
        line = ""
        line += (
            "<epsilon>'"
            + ",".join(
                [
                    ";".join(map(str, data["epsilon"]["epsilon"][:, i]))
                    for i in range(0, 3)
                ]
            )
            + "'</epsilon>"
        )
        line += (
            "<epsilon_rpa>"
            + ",".join(map(str, data["epsilon"]["epsilon_rpa"].flatten()))
            + "</epsilon_rpa>"
        )
        line += (
            "<epsilon_ion>"
            + ",".join(map(str, data["epsilon"]["epsilon_ion"].flatten()))
            + "</epsilon_ion>"
        )
        born_charges = data["born_charges"]
        atoms = vrun.all_structures[-1]
        spgg = Spacegroup3D(atoms)
        wycs = spgg._dataset["wyckoffs"]
        spg = str(spgg._dataset["number"])
        natoms = atoms.num_atoms
        combs = []
        line += "<born_effective_charge>"
        for k in range(natoms):
            comb = str(atoms.elements[k]) + "," + str(wycs[k])
            if comb not in combs:
                p = (
                    str(atoms.elements[k])
                    + ","
                    + str(wycs[k])
                    + ","
                    + ",".join(map(str, born_charges[k].flatten()))
                    + ";"
                )
                line += p
                combs.append(comb)
        line += "</born_effective_charge>"
        # print (line)
        vrun_eigs = data["phonon_eigenvalues"]
        pza = out.piezoelectric_tensor[1]
        pz = ",".join([";".join(map(str, pza[:, i])) for i in range(0, 6)])
        line += (
            '<dfpt_piezoelectric_tensor>"'
            + pz
            + '"</dfpt_piezoelectric_tensor>'
        )

        phonon_eigenvalues = out.phonon_eigenvalues
        phonon_eigenvectors = data["phonon_eigenvectors"]
        masses = data["masses"]
        born_charges = data["born_charges"]
        x, y = ir_intensity(
            phonon_eigenvectors=phonon_eigenvectors,
            phonon_eigenvalues=phonon_eigenvalues,
            masses=masses,
            born_charges=born_charges,
        )

        line += "<ir_intensity>'"
        line += ",".join(map(str, x)) + ";" + ",".join(map(str, y))
        line += "'</ir_intensity>"
        return line

    def main_lepsilon(self):
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-LEPSILON*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                data = self.dfpt_related(out=out, vrun=vrun)
            except Exception:
                print("Cannot find LEPSION data", folder)
                pass
        os.chdir(folder)
        info["main_lepsilon_info"] = data
        self.main_lepsilon_info = info
        return info

    def elastic_props(
        self, outcar="MAIN-ELASTIC-bulk@mp-149/OUTCAR",
    ):
        out = Outcar(outcar)
        cij = np.array(out.elastic_props()["cij"])
        d = ElasticTensor(cij).to_dict()
        line = ""
        line += (
            '<cij>"'
            + ";".join([",".join(map(str, cij[:, i])) for i in range(0, 6)])
            + '"</cij>'
        )
        # print("line", line)
        for i, j in d.items():
            line += "<" + str(i) + ">'" + str(j) + "'</" + str(i) + ">"
        totdos = outcar.replace("OUTCAR", "total_dos.dat")
        if os.path.isfile(totdos):

            mesh_yaml = outcar.replace("OUTCAR", "mesh.yaml")
            band_yaml = outcar.replace("OUTCAR", "band.yaml")
            with open(mesh_yaml, "r") as f:
                doc = yaml.load(f)
            nmodes = doc["phonon"][0]["band"]
            ph_modes = []

            for p in nmodes:
                ph_modes.append(p["frequency"])
            ph_modes = sorted(set(ph_modes))
            f = open(totdos, "r")
            freq = []
            pdos = []
            for lines in f.readlines():
                if not str(lines.split()[0]).startswith("#"):
                    #   print (lines)
                    # else:
                    freq.append(float(lines.split()[0]))
                    pdos.append(float(lines.split()[1]))
            line += (
                "<phonon_modes>"
                + ",".join(map(str, [round(i, 2) for i in ph_modes]))
                + "</phonon_modes>"
            )
            line += (
                "<phonon_dos_frequencies>'"
                + ",".join(map(str, freq))
                + "'</phonon_dos_frequencies>"
            )
            line += (
                "<phonon_dos_intensity>'"
                + ",".join(map(str, pdos))
                + "'</phonon_dos_intensity>"
            )
            from jarvis.io.phonopy.outputs import bandstructure_plot

            frequencies, distances, labels, label_points = bandstructure_plot(
                band_yaml
            )

            tmp = ""
            for i in range(np.array(frequencies).shape[1]):
                tmp += ",".join(map(str, np.array(frequencies)[:, i])) + ";"
            line += (
                "<phonon_bandstructure_distances>'"
                + ",".join(map(str, distances))
                + "'</phonon_bandstructure_distances>"
            )
            line += (
                "<phonon_bandstructure_frequencies>'"
                + tmp
                + "'</phonon_bandstructure_frequencies>"
            )
            line += (
                "<phonon_bandstructure_labels>'"
                + ",".join(map(str, labels))
                + "'</phonon_bandstructure_labels>"
            )
            line += (
                "<phonon_bandstructure_label_points>'"
                + ",".join(map(str, label_points))
                + "'</phonon_bandstructure_label_points>"
            )
        return line

    def main_elastic(self):
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-ELAST*.json"):
            try:
                folder = i.split(".json")[0]
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                data = self.elastic_props(out)
            except Exception:
                print("Cannot get ELAST folder", folder)
                pass
        os.chdir(folder)
        info["main_elastic_info"] = data
        self.main_elastic_info = info
        return info

    def boltztrap_data(
        self,
        path="MAIN-RELAX-bulk@mp-149/boltztrap",
        temperature=600,
        doping=1e20,
    ):
        info = {}

        all_data = BoltzTrapOutput(path).to_dict()
        # info["all_data"] = all_data
        small_p = all_data["condtens_fixdoping"]["p"][temperature]
        small_n = all_data["condtens_fixdoping"]["n"][temperature]
        for i, j in small_p.items():
            if j["N_cm3"] == doping:
                tmp = j
                pseeb = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(
                            tmp["seeb"].reshape(3, 3) * 1e6
                        )
                    ]
                )
                pcond = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(tmp["cond"].reshape(3, 3))
                        / 1e14
                    ]
                )
                ppf = pseeb ** 2 * pcond / 1e6
                pkappa = np.linalg.eigvals(tmp["kappa"].reshape(3, 3))

                info["pseeb"] = pseeb.tolist()
                info["pcond"] = pcond.tolist()
                info["ppf"] = ppf.tolist()
                info["pkappa"] = pkappa.tolist()

        for i, j in small_n.items():
            if j["N_cm3"] == -1 * doping:
                tmp = j
                nseeb = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(
                            tmp["seeb"].reshape(3, 3) * 1e6
                        )
                    ]
                )
                ncond = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(tmp["cond"].reshape(3, 3))
                        / 1e14
                    ]
                )
                npf = nseeb ** 2 * ncond / 1e6
                nkappa = np.linalg.eigvals(tmp["kappa"].reshape(3, 3))

                info["nseeb"] = np.array([np.real(r) for r in nseeb.tolist()])
                info["ncond"] = np.array([np.real(r) for r in ncond.tolist()])
                info["npf"] = npf.tolist()
                info["nkappa"] = nkappa.tolist()
        line = ""
        for i, j in info.items():
            if isinstance(j, list):
                tmp = ";".join(map(str, j))
            else:
                tmp = j
            line += (
                "<"
                + str(i)
                + ">'"
                + ",".join(map(str, j))
                + "'</"
                + str(i)
                + ">"
            )
        # print (line)
        # print (all_data)
        return line

    def main_boltz_data(self):
        folder = self.folder
        info = {}
        data = ""
        os.chdir(folder)
        for i in glob.glob("MAIN-RELAX*.json"):
            try:
                folder = i.split(".json")[0]
                path = os.getcwd() + "/" + folder + "/boltztrap"
                data = self.boltztrap_data(path=path)
            except Exception:
                print("Cannot find boltztrap folder", folder)
                pass
        os.chdir(folder)
        info["boltztrap_info"] = data
        self.boltztrap_info_dict = info
        return info

    def image_to_string(
        self,
        img_path="/cluster/users/knc6/justback/2DSTM/PNG_JSON3/JVASP-60776_mp-19795_pos.jpg",
    ):
        # 2D array only
        fig = imread(img_path)
        line = (
            str(fig.shape[0])
            + "_"
            + str(fig.shape[1])
            + "_"
            + ",".join(map(str, fig[:, :, 0].flatten()))
        )
        return line

    def main_stm(self):
        folder = self.folder
        # for i in glob.glob()
        pass

    def efg_tensor(self):
        main_folder = self.folder
        info = {}
        efg_dat = ""
        os.chdir(main_folder)
        for i in glob.glob("MAIN-LEFG*.json"):
            try:
                folder = i.split(".json")[0]
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                atoms = Atoms.from_poscar(out.replace("OUTCAR", "POSCAR"))
                out = Outcar(out)
                spgg = Spacegroup3D(atoms)
                wycs = spgg._dataset["wyckoffs"]
                natoms = atoms.num_atoms
                combs = []
                efg = out.efg_raw_tensor
                efg_dat = ""
                for k in range(natoms):
                    comb = str(atoms.elements[k]) + "," + str(wycs[k])
                    if comb not in combs:
                        efg_dat = (
                            efg_dat
                            + comb
                            + ","
                            + ",".join(map(str, efg[k].flatten()))
                            + ";"
                        )
                        combs.append(comb)
            except Exception:
                print("Cannot get LEFG folder", folder)
                pass

        info["efg_raw_tensor"] = efg_dat
        self.efg_raw_tensor_info = info
        os.chdir(main_folder)
        return info

    def write_xml(self, filename="temp.xml"):

        with open(filename, "w") as f:
            line = '<?xml version="1.0" encoding="UTF-8"?>\n<?xml-stylesheet type="text/xsl" href="temp-1002.xsl"?>\n<basic_info>'
            f.write(str(line))
            line = (
                "<convergence_info>"
                + stringdict_to_xml(self.encut_kp())
                + "</convergence_info>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_relax_info>"
                + stringdict_to_xml(self.basic_info())
                + "</main_relax_info>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_band>"
                + stringdict_to_xml(self.main_band(), enforce_string=False)
                + "</main_band>"
            )
            f.write(str(line))
            line = (
                "<main_optics_semilocal>"
                + stringdict_to_xml(self.main_optics_semilocal())
                + "</main_optics_semilocal>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_optics_mbj>"
                + stringdict_to_xml(self.main_optics_mbj())
                + "</main_optics_mbj>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_elastic>"
                + stringdict_to_xml(self.main_elastic())
                + "</main_elastic>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_boltz>"
                + stringdict_to_xml(self.main_boltz_data())
                + "</main_boltz>"
                + "\n"
            )
            f.write(str(line))
            line = stringdict_to_xml(self.main_lepsilon()) + "\n"
            f.write(str(line))
            line = stringdict_to_xml(self.main_soc_spillage()) + "\n"
            f.write(str(line))
            line = (
                stringdict_to_xml(self.efg_tensor(), enforce_string=True)
                + "\n"
            )
            f.write(str(line))
            f.write("</basic_info>\n")
            """
                line = "<stm_image>" + image_to_string() + "</stm_image>"+ "\n"
                f.write(str(line))
                line ='</basic_info>'
                f.write(str(line))
                """


VaspToApiXmlSchema(
    folder="/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO"
).write_xml(filename="JVASP-1002.xml")

# directories = ["/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO"]
# filenames = ["JVASP-1002.xml"]
# main_compile(directories=directories, filenames=filenames)
