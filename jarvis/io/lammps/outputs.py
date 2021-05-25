"""Function to analze LAMMPS output."""
import numpy as np
import glob
import os
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.io.lammps.inputs import LammpsData
from jarvis.io.phonopy.outputs import bandstructure_plot, total_dos


def parse_potential_mod(mod="potential.mod"):
    """Parse potentials.mod input file."""
    info = {}
    f = open(mod, "r")
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        if "pair_style" in i:
            pair_style = i.split("pair_style")[1]
        if "pair_coeff" in i:
            pair_coeff = i.split()[3].split("/")[-1]
            elements = i.split()[4:]

    info["pair_style"] = pair_style
    info["pair_coeff"] = pair_coeff
    info["elements"] = elements
    return info


def read_data(data=None, ff=None, element_order=[]):
    """
    Read LAMMPS data file.

    Args:
        data: data file path

        ff: potential.mod/potential information file path

    Returns:
          Atoms object
    """
    lmp_data = LammpsData()
    return lmp_data.read_data(
        filename=data,
        element_order=element_order,
        potential_file=ff,
        verbose=False,
    )


def parse_log(log="log.lammps"):
    """
    Analyzes log.lammps file.

    Please note, the output format heavily depends on the input file
    A generic input is taken here.
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
    Et = ""
    info = {}
    try:
        logfile = open(log, "r")
        lines = logfile.read().splitlines()
        for i, line in enumerate(lines):
            if "Loop time of" in line:
                # toten = float(lines[i - 1].split()[12])
                press = float(lines[i - 1].split()[2])
                press = float(press) * 0.0001
                denom = float(lines[i - 1].split()[17])
                en = float(lines[i - 1].split()[12]) / denom
                break
        logfile.close()
    except Exception:
        print("Cannot parse energy and pressure information.", log)
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
            if 'print "Elastic Constant C15all = ${C15all} ${cunits}"' in line:
                c15 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
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
        cij = np.array(
            [
                [c11, c12, c13, c14, c15, c16],
                [c12, c22, c23, c24, c25, c26],
                [c13, c23, c33, c34, c35, c36],
                [c14, c24, c34, c44, c45, c46],
                [c15, c25, c35, c45, c55, c56],
                [c16, c26, c36, c46, c56, c66],
            ],
            dtype="float",
        )
        Et = ElasticTensor(et_tensor=cij).to_dict()
    except Exception:
        print("Cannot parse log.lammps", log)
        pass
    info["energy_per_atom"] = en
    info["system_pressure"] = press
    info["elastic_tensor"] = Et
    return info


def get_chem_pot(all_json_data={}):
    """Get minumum energy for each elemental specie."""
    all_possible_species = []
    elemental_energies = {}
    for i, j in all_json_data.items():
        #  print (i)
        if (
            i != "jid"
            and i != "source_folder"
            and i != "bulk_data"
            and i != "bulk_energy_per_atom"
        ):
            for el in j["final_str"].elements:
                if el not in all_possible_species:
                    all_possible_species.append(el)
            if len(j["final_str"].uniq_species) == 1:
                tmp_element = j["final_str"].uniq_species[0]
                tmp_energy = j["energy_per_atom"]
                elemental_energies.setdefault(tmp_element, []).append(
                    [tmp_energy, i]
                )
    chem_pot = {}
    for i, j in elemental_energies.items():
        chem_pot[i] = sorted(j, key=lambda x: x[0])[0][0]
    if len(all_possible_species) != len(elemental_energies.keys()):
        raise ValueError(
            "Error: Chemical potential for all elements are not available"
        )
        raise ValueError(all_possible_species, elemental_energies.keys())
    return chem_pot


def parse_folder(path="bulk@mp-1487_fold/bulk@mp-1487",):
    """Parse individual LAMMPS run."""
    info = {}
    ff = os.path.join(path, "potential.mod")
    pot = parse_potential_mod(ff)
    initial_str = read_data(data=os.path.join(path, "data"), ff=ff)
    final_str = read_data(data=os.path.join(path, "data0"), ff=ff)
    log_path = os.path.join(path, "log.lammps")
    pair_style = pot["pair_style"]
    pair_coeff = pot["pair_coeff"]
    info["pair_style"] = pair_style
    info["pair_coeff"] = pair_coeff
    info["initial_str"] = initial_str
    info["final_str"] = final_str

    dat = parse_log(log_path)
    info["energy_per_atom"] = dat["energy_per_atom"]
    info["system_pressure"] = dat["system_pressure"]
    info["elastic_tensor"] = dat["elastic_tensor"]

    return info


def parse_material_calculation_folder(path="bulk@mp-1487_fold", jid="x"):
    """
    Parse individual LAMMPS material run.

    with optimization, vacancy, phonon, surface etc.
    """
    cwd = os.getcwd()
    jid_file = os.path.join(path, "JARVISFF-ID")
    if os.path.exists(jid_file):
        f = open(jid_file, "r")
        lines = f.read().splitlines()
        f.close()
        jid = lines[0]
    print(path)
    info = {}
    info["jid"] = jid
    info["source_folder"] = path
    bulk_energy_per_atom = ""
    for i in glob.glob(path + "/*.json"):
        try:
            json_file_name = i.split("/")[-1]
            json_file_path = i.split(".json")[0]
            print("json_file_name", json_file_name)
            print("json_file_path", json_file_path)
            fold_path = os.path.join(path, json_file_path)
            tmp_info = parse_folder(fold_path)
            info[json_file_name] = tmp_info
            if (
                "bulk" in json_file_name
                and "cellmax" not in json_file_name
                and "sbulk" not in json_file_name
            ):
                bulk_energy_per_atom = tmp_info["energy_per_atom"]
                info["bulk_data"] = tmp_info
        except Exception:
            print("Error cannot parse bulk file.")
            pass
    info["bulk_energy_per_atom"] = bulk_energy_per_atom
    # Each element has energy_per_atom,
    # system_pressure,elastic_tensor,final_str,initial_str
    print("Found", len(info.keys()), "folders")
    try:
        chem_pot = get_chem_pot(info)
        info["chem_pot"] = chem_pot
    except Exception:
        print("Seems like didnt run vacancy calcs.")
        pass
    # print (chem_pot)
    vacancy_info = []
    surface_info = []
    try:
        for i, j in info.items():
            if "vacancy" in i:
                try:
                    element_vacant = i.split("-")[-1].split("@")[0]
                    perfect_energy = (
                        j["final_str"].num_atoms + 1
                    ) * bulk_energy_per_atom
                    defect_energy = (
                        j["final_str"].num_atoms * j["energy_per_atom"]
                    )
                    mu = chem_pot[element_vacant]
                    Ef = defect_energy - perfect_energy + mu
                    mult = i.split("mult-")[1].split("_")[0]
                    vacancy_info.append([element_vacant, mult, Ef])
                    # print ('Ef',i,element_vacant, Ef,mult)
                except Exception:
                    print("Error parsing vacancy.", i)
                    pass
            if "Surf" in i:
                try:
                    perfect_energy = (
                        j["final_str"].num_atoms
                    ) * bulk_energy_per_atom
                    defect_energy = (
                        j["final_str"].num_atoms * j["energy_per_atom"]
                    )
                    m = j["final_str"].lattice.matrix
                    area = np.linalg.norm(np.cross(m[0], m[1]))
                    Ef = 16 * (defect_energy - perfect_energy) / (2 * area)
                    surf_name = i.split("@")[0].split("-")[1].replace("_", " ")
                    surface_info.append([surf_name, Ef])
                except Exception:
                    print("Error parsing surface.", i)
                    pass
                # print (i,Ef)
    except Exception:
        print("Error parsing vacancy-surface.")
        pass
    # print ('JID',jid)
    info["surface_info"] = surface_info
    info["vacancy_info"] = vacancy_info
    dos_file = os.path.join(path, "Phonon", "total_dos.dat")
    # mesh_file = os.path.join(path, "Phonon", "mesh.yaml")
    band_file = os.path.join(path, "Phonon", "band.yaml")
    if os.path.exists(dos_file):
        dos_freq, dos_intensity = total_dos(tot_dos=dos_file)
        info["dos_freq"] = dos_freq
        info["dos_intensity"] = dos_intensity
    if os.path.exists(band_file):
        (
            band_frequencies,
            band_distances,
            band_labels,
            band_label_points,
        ) = bandstructure_plot(band_file)
        info["band_frequencies"] = band_frequencies
        info["band_distances"] = band_distances
        info["band_labels"] = band_labels
        info["band_label_points"] = band_label_points
    os.chdir(cwd)
    return info


def parse_full_ff_folder(path="Mishin-Ni-Al-2009.eam.alloy_nist",):
    """Parse complete FF calculation folder."""
    cwd = os.getcwd()
    os.chdir(path)
    print("path", path)
    tmp_path = path + "/*_fold"
    for i in glob.glob(tmp_path):
        print(i)
        tmp_fold = os.path.join(path, i)
        info = parse_material_calculation_folder(tmp_fold)
    os.chdir(cwd)
    return info


def analyze_log(log="log.lammps"):
    """
    Analyzes log.lammps file.

    Please note, the output format heavily depends on the input file
    A generic input is taken here.
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
    # c15 = 0
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
                denom = float(lines[i - 1].split()[17])
                en = float(lines[i - 1].split()[12]) / denom
                break
        logfile.close()
    except Exception:
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
    except Exception:
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
        # round(float(c15), 1),
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


def read_dump(data=None):
    """Read LAMMPS dump file."""
    f = open(data, "r")
    lines = f.read().splitlines()
    for i, line in enumerate(lines):
        if "NUMBER OF ATOMS" in line:
            natoms = int(lines[i + 1].split()[0])
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    z = np.zeros((natoms))
    coords = list()  # np.zeros((natoms))
    for i, line in enumerate(lines):
        if "ITEM: ATOMS" in line:
            for j in range(0, natoms):
                x[j] = (lines[i + j + 1]).split()[1]
                y[j] = (lines[i + j + 1]).split()[2]
                z[j] = (lines[i + j + 1]).split()[3]
                coords.append([x[j], y[j], z[j]])
    f.close()
    prop = np.asarray(coords)
    return prop


# p=read_data()
# print (p)
# parse_potential_mod()
"""
if __name__ == "__main__":
    lg = "log.lammps"
    x = analyze_log(lg)
    print(x)
"""
