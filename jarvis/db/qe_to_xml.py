"""Module to generate XML file for QE calculation."""
import gzip
import os
import xmltodict
from jarvis.core.atoms import get_supercell_dims
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.utils import stringdict_to_xml, xml_to_dict, array_to_string
from jarvis.io.qe.outputs import ProjHamXml, DataFileSchema
from jarvis.analysis.thermodynamics.energetics import (
    form_enp,
    get_unary_qe_tb_energy,
)

# import pprint


def parse_material_calculation_folder(
    path="Ag_Al_kspace/dimer.in_vnscf_coords_min_2",
    source="JARVIS-DFT-QE",
):
    """Parse QE calculation folder."""
    save_path = os.path.join(path, "qe.save")
    info = {}
    if not os.path.exists(save_path):
        raise ValueError("Doesnt exist:", save_path)
    jid_path = os.path.join(path, "JARVIS-ID")
    keys = [
        "jid",
        "source_folder",
        "final_energy",
        "final_energy_breakdown",
        # "initial_structure",
        "natoms",
        "final_spacegroup_number",
        "final_spacegroup_symbol",
        "final_pointgroup_symbol",
        "final_crystal_system",
        # "initial_spacegroup_number",
        # "initial_spacegroup_symbol",
        # "initial_pointgroup_symbol",
        # "initial_crystal_system",
        "energy_per_atom",
        "elements",
        "formula",
        "number_uniq_species",
        "final_structure",
        "efermi",
        "functional",
        "indir_gap",
        "nkpts",
        "qe_version",
        "is_spin_polarized",
        "is_spin_orbit",
        "data_source",
        "forces",
        "stress",
        "nelec",
        "f_enp",
        "dos",
    ]
    for i in keys:
        info[i] = "na"
    if os.path.exists(jid_path):
        f = open(jid_path, "r")
        jid = f.read().splitlines()[0]
        f.close()

        info["jid"] = jid
        info["tmp_id"] = "'" + jid + "'"
        info["source_folder"] = "'" + path + "'"

    for i in os.listdir(save_path):
        if "data-file-schema" in i:
            if ".gz" in i:
                dat_path = os.path.join(save_path, i)
                f = gzip.open(dat_path, "rb")
                file_content = f.read()
                data = xmltodict.parse(file_content)
                if "output" in data["qes:espresso"]:
                    set_key = "output"
                elif "step" in data["qes:espresso"]:
                    set_key = "step"
                else:
                    raise ValueError("Inconsisten QE version.")
                data_schm = DataFileSchema(data=data, set_key=set_key)
            else:
                dat_path = os.path.join(save_path, i)
                data = xml_to_dict(dat_path)
                if "output" in data["qes:espresso"]:
                    set_key = "output"
                elif "step" in data["qes:espresso"]:
                    set_key = "step"
                else:
                    raise ValueError("Inconsisten QE version.")
                data_schm = DataFileSchema(data=data, set_key=set_key)
            info["final_energy"] = data_schm.final_energy
            for ii, jj in data_schm.final_energy_breakdown.items():
                info[ii] = jj
            final_strt = data_schm.final_structure
            info["final_structure"] = (
                "'"
                + data_schm.final_structure.get_string().replace("\n", "\\n")
                + "'"
            )
            try:
                print("data_schm.final_energy", data_schm.final_energy)
                f_enp = form_enp(
                    atoms=final_strt,
                    chem_pots=get_unary_qe_tb_energy(),
                    total_energy=data_schm.final_energy,
                )
                info["f_enp"] = round(f_enp, 4)
            except Exception as exp:
                print(exp)
                pass
            dim = get_supercell_dims(final_strt)
            info["xyz"] = (
                '"'
                + str(
                    final_strt.make_supercell_matrix(dim).get_xyz_string
                ).replace("\n", "\\n")
                + '"'
            )
            info["natoms"] = final_strt.num_atoms
            final_spg = Spacegroup3D(final_strt)
            info["final_spacegroup_number"] = final_spg.space_group_number
            info["final_spacegroup_symbol"] = final_spg.space_group_symbol
            info["final_pointgroup_symbol"] = final_spg.point_group_symbol
            info["final_crystal_system"] = final_spg.crystal_system
            final_lat_params = final_strt.lattice.parameters
            info["final_a"] = round(final_lat_params[0], 2)
            info["final_b"] = round(final_lat_params[1], 2)
            info["final_c"] = round(final_lat_params[2], 2)
            info["final_alpha"] = round(final_lat_params[3], 2)
            info["final_beta"] = round(final_lat_params[4], 2)
            info["final_gamma"] = round(final_lat_params[5], 2)
            info["final_density"] = round(final_strt.density, 3)

            initial_strt = data_schm.initial_structure
            initial_lat_params = initial_strt.lattice.parameters
            info["initial_a"] = round(initial_lat_params[0], 2)
            info["initial_b"] = round(initial_lat_params[1], 2)
            info["initial_c"] = round(initial_lat_params[2], 2)
            info["initial_alpha"] = round(initial_lat_params[3], 2)
            info["initial_beta"] = round(initial_lat_params[4], 2)
            info["initial_gamma"] = round(initial_lat_params[5], 2)
            info["initial_density"] = round(initial_strt.density, 3)
            """
            initial_spg = Spacegroup3D(initial_strt)
            info["initial_spacegroup_number"] = initial_spg.space_group_number
            info["initial_spacegroup_symbol"] = initial_spg.space_group_symbol
            info["initial_pointgroup_symbol"] = initial_spg.point_group_symbol
            info["initial_crystal_system"] = initial_spg.crystal_system
            info["initial_structure"] = (
                "'"
                + data_schm.initial_structure.get_string().replace("\n", "\\n")
                + "'"
            )
            """

            info["energy_per_atom"] = round(
                float(data_schm.energy_per_atom), 4
            )

            info["cif"] = ""
            try:
                cif_filename = os.path.join(os.getcwd(), "atoms.cif")
                final_strt.write_cif(filename="atoms.cif")
                cif_f = open(cif_filename, "r")
                cif_lines = cif_f.read()  # .splitlines()
                cif_f.close()
                os.remove(cif_filename)
                info["cif"] = '"' + str(cif_lines).replace("\n", "\\n") + '"'
            except Exception as exp:
                print("CIF exception", exp)
                pass
            info["elements"] = ",".join(final_strt.uniq_species)
            info["tmp_elements"] = (
                "'" + ",".join(final_strt.uniq_species) + "'"
            )
            info["formula"] = final_strt.composition.reduced_formula
            info["tmp_formula"] = (
                "'" + final_strt.composition.reduced_formula + "'"
            )
            info["number_uniq_species"] = len(final_strt.uniq_species)
            info["efermi"] = round(data_schm.efermi, 4)
            info["functional"] = data_schm.functional
            info["indir_gap"] = round(data_schm.indir_gap, 3)
            info["nkpts"] = data_schm.nkpts
            info["nelec"] = data_schm.nelec
            info["qe_version"] = data_schm.qe_version
            info["is_spin_polarized"] = data_schm.is_spin_polarized
            info["is_spin_orbit"] = data_schm.is_spin_orbit
            info["data_source"] = source
            info["forces"] = (
                "'"
                + array_to_string(
                    [" ".join(map(str, j)) for j in data_schm.forces]
                )
                + "'"
            )
            info["stress"] = (
                "'"
                + array_to_string(
                    [" ".join(map(str, j)) for j in data_schm.stress]
                )
                + "'"
            )
    try:
        projham_xml_path_gz = os.path.join(path, "projham_K.xml.gz")
        projham_xml_path = os.path.join(path, "projham_K.xml")
        if os.path.exists(projham_xml_path):
            energies, dos, pdos, names = ProjHamXml(
                filename=projham_xml_path
            ).dos()
            line = (
                "<edos_energies>'"
                + ",".join(map(str, energies))
                + "'</edos_energies>"
            )
            line += (
                "<total_edos_up>'"
                + ",".join(map(str, dos))
                + "'</total_edos_up>"
            )
            line += "<elemental_dos>"
            for n, nn in enumerate(names):
                pdos_n = pdos[:, n]
                line += str(nn) + "_" + ",".join(map(str, pdos_n)) + ";"
            line += "</elemental_dos>"

            info["dos"] = line
        elif os.path.exists(projham_xml_path_gz):
            energies, dos, pdos, names = ProjHamXml(
                filename=projham_xml_path_gz
            ).dos()
            line = (
                "<edos_energies>'"
                + ",".join(map(str, energies))
                + "'</edos_energies>"
            )
            line += (
                "<total_edos_up>'"
                + ",".join(map(str, dos))
                + "'</total_edos_up>"
            )
            line += "<elemental_dos>'"
            for n, nn in enumerate(names):
                pdos_n = pdos[:, n]
                line += str(nn) + "_" + ",".join(map(str, pdos_n)) + ";"
            line += "'</elemental_dos>"
            info["dos"] = line
        else:
            info["dos"] = "na"
    # print(pprint.pprint(info))
    except Exception as exp:
        print(exp)
        pass
    return info


def write_xml(path="", filename="temp.xml"):
    """Write XML file."""
    basic = parse_material_calculation_folder(path=path)
    f = open(filename, "w")
    basic_xml = stringdict_to_xml(basic)
    line = ""
    line = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line += '<?xml-stylesheet type="text/xsl" '
    line += 'href="jarvisqe.xsl"?>\n'

    line += "<basic_info>"

    line += basic_xml
    line += "</basic_info>"
    f.write(line)
    f.close()


"""
#path='/rk2/knc6/UniveralTB/julia_data/TEST_ELS/POSCAR_JVASP-1002'
#path='/rk2/knc6/UniveralTB/julia_data/TEST_BIN/POSCAR_JVASP-10016'
if __name__ == "__main__":
    #path = "/rk2/knc6/UniveralTB/julia_data/
    #Si_C_kspace/POSCAR_tio2_rutile_vnscf_vol_1"
    write_xml(path=path)
"""
