"""Module to generate XML file for LAMMPS calculation."""
import numpy as np
from jarvis.core.atoms import get_supercell_dims
from jarvis.db.jsonutils import loadjson
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.utils import stringdict_to_xml

mp_jv = loadjson("/rk2/knc6/DB/MP/mp_jv_id.json")


def get_jvid(mp=""):
    """Get JARVIS-ID for MPID."""
    jvid = ""
    try:
        jvid = "'" + ",".join(mp_jv[mp]) + "'"
    except Exception as exp:
        print("No JID", exp)
        pass
    return jvid


def basic_data(data={}, source="JARVIS-FF-LAMMPS"):
    """Get basic data for table."""
    info = {}
    info["id"] = data["jid"]
    info["source_folder"] = data["source_folder"]
    info["tmp_source_folder"] = "'" + data["source_folder"] + "'"
    info["tmp_id"] = "'" + data["jid"] + "'"
    ref = data["source_folder"].split("/")[-1].split("@")[1].split("_")[0]
    info["ref"] = "'" + ref + "'"
    info["jvid"] = get_jvid(ref)
    final_atoms = data["bulk_data"]["final_str"]
    initial_atoms = data["bulk_data"]["initial_str"]
    info["formula"] = final_atoms.composition.reduced_formula
    info["tmp_formula"] = "'" + final_atoms.composition.reduced_formula + "'"
    info["elements"] = ",".join(final_atoms.uniq_species)
    info["tmp_elements"] = "'" + ",".join(final_atoms.uniq_species) + "'"
    info["number_uniq_species"] = len(final_atoms.uniq_species)
    info["data_source"] = source
    info["tmp_data_source"] = "'" + source + "'"
    info["pair_style"] = data["bulk_data"]["pair_style"]
    info["pair_coeff"] = data["bulk_data"]["pair_coeff"]
    # info["pair_coeff"] = (
    #    '<a href="http://www.ctcms.nist.gov/~knc6/DOWNLOADS/'
    #    + data["bulk_data"]["pair_coeff"]
    #    + '.zip">'
    #    + data["bulk_data"]["pair_coeff"]
    #    + "</a>"
    # )
    info["energy_per_atom"] = round(
        float(data["bulk_data"]["energy_per_atom"]), 3
    )
    info["pressure"] = round(float(data["bulk_data"]["system_pressure"]), 3)

    initial_spg = Spacegroup3D(initial_atoms)
    final_spg = Spacegroup3D(final_atoms)
    info["initial_spacegroup_number"] = initial_spg.space_group_number
    info["initial_spacegroup_symbol"] = initial_spg.space_group_symbol
    info["initial_pointgroup_symbol"] = initial_spg.point_group_symbol
    info["initial_crystal_system"] = initial_spg.crystal_system
    initial_lat_params = initial_atoms.lattice.parameters
    info["initial_a"] = round(initial_lat_params[0], 2)
    info["initial_b"] = round(initial_lat_params[1], 2)
    info["initial_c"] = round(initial_lat_params[2], 2)
    info["initial_alpha"] = round(initial_lat_params[3], 2)
    info["initial_beta"] = round(initial_lat_params[4], 2)
    info["initial_gamma"] = round(initial_lat_params[5], 2)
    info["initial_density"] = round(initial_atoms.density, 3)
    dim = get_supercell_dims(final_atoms)
    info["xyz"] = (
        '"'
        + str(final_atoms.make_supercell_matrix(dim).get_xyz_string).replace(
            "\n", "\\n"
        )
        + '"'
    )
    info["final_str"] = (
        '"' + str(final_atoms.get_string()).replace("\n", "\\n") + '"'
    )
    info["initial_str"] = (
        '"' + str(initial_atoms.get_string()).replace("\n", "\\n") + '"'
    )
    final_lat_params = final_atoms.lattice.parameters
    info["final_a"] = round(final_lat_params[0], 2)
    info["final_b"] = round(final_lat_params[1], 2)
    info["final_c"] = round(final_lat_params[2], 2)
    info["final_alpha"] = round(final_lat_params[3], 2)
    info["final_beta"] = round(final_lat_params[4], 2)
    info["final_gamma"] = round(final_lat_params[5], 2)
    info["final_density"] = round(final_atoms.density, 3)

    info["final_spacegroup_number"] = final_spg.space_group_number
    info["final_spacegroup_symbol"] = final_spg.space_group_symbol
    info["final_pointgroup_symbol"] = final_spg.point_group_symbol
    info["final_crystal_system"] = final_spg.crystal_system

    et = ""
    # print(data["bulk_data"]["elastic_tensor"]["raw_et_tensor"])
    try:
        if data["bulk_data"]["elastic_tensor"] != "":
            cij = np.round(
                (data["bulk_data"]["elastic_tensor"]["raw_et_tensor"]), 2
            )
            et = (
                '<cij>"'
                + ";".join(
                    [",".join(map(str, cij[:, i])) for i in range(0, 6)]
                )
                + '"</cij>'
            )

    except Exception as exp:
        print("Cannot obtain elastic tensor data.", info["source_folder"], exp)
        pass

    info["elastic_tensor"] = et
    vacancy_line = ""
    if len(data["vacancy_info"]) != 0:
        for i in data["vacancy_info"]:
            vacancy_line += i[0] + "," + i[1] + "," + str(round(i[2], 3)) + ";"
    info["vacancy_info"] = '"' + vacancy_line + '"'

    surf_line = ""
    if len(data["surface_info"]) != 0:
        for i in data["surface_info"]:
            surf_line += i[0] + "," + str(round(i[1], 3)) + ";"
    info["surface_info"] = '"' + surf_line + '"'

    phonon_band_line = ""
    try:
        # Band
        frequencies = data["band_frequencies"]
        distances = data["band_distances"]
        labels = data["band_labels"]
        label_points = data["band_label_points"]

        tmp = ""
        for i in range(np.array(frequencies).shape[1]):
            tmp += ",".join(map(str, np.array(frequencies)[:, i])) + ";"

        phonon_band_line += (
            "<phonon_bandstructure_distances>'"
            + ",".join(map(str, distances))
            + "'</phonon_bandstructure_distances>"
        )

        phonon_band_line += (
            "<phonon_bandstructure_frequencies>'"
            + tmp
            + "'</phonon_bandstructure_frequencies>"
        )
        phonon_band_line += (
            "<phonon_bandstructure_labels>'"
            + ",".join(map(str, labels))
            + "'</phonon_bandstructure_labels>"
        )
        phonon_band_line += (
            "<phonon_bandstructure_label_points>'"
            + ",".join(map(str, label_points))
            + "'</phonon_bandstructure_label_points>"
        )
    except Exception as exp:
        print("Cannot obtain phonon band data.", exp)

    # Comment until 4 MB text size error
    info["phonon_band_line"] = phonon_band_line
    phonon_dos_line = ""
    try:
        freq = data["dos_freq"]
        pdos = data["dos_intensity"]
        phonon_dos_line += (
            "<phonon_dos_frequencies>'"
            + ",".join(map(str, freq))
            + "'</phonon_dos_frequencies>"
        )

        phonon_dos_line += (
            "<phonon_dos_intensity>'"
            + ",".join(map(str, pdos))
            + "'</phonon_dos_intensity>"
        )

    except Exception:
        print("Cannot obtain phonon dod data.")
        pass
    info["phonon_dos_line"] = phonon_dos_line

    return info


def write_xml(data={}, filename="temp.xml"):
    """Write XML file."""
    basic = basic_data(data)
    # name = data["jid"] + ".xml"
    f = open(filename, "w")
    basic_xml = stringdict_to_xml(basic)
    line = ""
    line = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line += '<?xml-stylesheet type="text/xsl" '
    line += 'href="jarvisff.xsl"?>\n'

    line += "<basic_info>"

    line += basic_xml
    line += "</basic_info>"
    f.write(line)
    f.close()


"""
if __name__ == "__main__":
    from jarvis.io.lammps.outputs import parse_material_calculation_folder

    fold = '/rk2/knc6/JARVIS-FF/ALLOY8/
    Mishin-Ni-Al-2009.eam.alloy_nist/bulk@mp-23_fold'
    x = parse_material_calculation_folder(fold)
    write_xml(x)
    print('src',x["source_folder"].split('/')[-1].split('@')[1].split('_')[0])
"""
