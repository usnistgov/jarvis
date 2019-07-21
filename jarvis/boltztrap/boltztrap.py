"""
Run and analyze BoltzTrap
"""
from __future__ import unicode_literals, print_function
from shutil import which
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer, BoltztrapRunner
from pymatgen.electronic_structure.plotter import BSPlotter, BoltztrapPlotter
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import json, os
from monty.serialization import MontyEncoder, loadfn, MontyDecoder


def boltz_run(vrun=""):
    """
  Helper function to run and store Bolztrap 
  Please note high denske k-point mesh is needed
  The automatic k-point convergene in JARVIS-DFT is generally good enough
  
  Args:
      vrun: path to vasprun.xml
  Returns:
        BoltztrapAnalyzer object
  """
    if which("x_trans") is None or which("BoltzTraP") is None:
        print("Please install BoltzTrap")
    v = Vasprun(vrun)
    kp = vrun.replace("vasprun.xml", "KPOINTS")
    out = vrun.replace("vasprun.xml", "OUTCAR")
    for line in open(out, "r"):
        if "NELECT" in line:
            nelect = float(line.split()[2])
    bs = v.get_band_structure(kp, line_mode=False)
    brun = BoltztrapRunner(bs, nelect)
    path = vrun.split("vasprun.xml")[0]
    folder = str(path) + str("/boltztrap")
    out_trans = str(folder) + str("/boltztrap.outputtrans")
    if not os.path.exists(out_trans):
        brun.run(path_dir=path)
        print("doesnt exist", out_trans)
    bana = BoltztrapAnalyzer.from_files(folder)
    return bana


def get_prop(bzana, prop="seebeck", dop_type="n", temp=700, dop=1e20, output="eigs"):
    """
    Helper function to obtain thermoelectric properties at specific values of temp,doping
    Args:
       bzana: BoltztrapAnalyzer object
       prop: one of 'seebeck','power_factor', 'electronic_conductivity', 'el_thermal_conductivity', 'zt'
       dop_type: 'n' or 'p' 
       temp: temperature in K
       dop: doping conc. in /cm3
       output: 'eigs' or 'tensor'
    Returns:
          val: desired values
    """
    doping = bzana.doping["p"]
    index = None
    for i, ii in enumerate(doping):
        if ii == dop:
            index = i
    if index is not None:
        if prop == "seebeck":
            val = bzana.get_seebeck(output=output)[dop_type][temp][index]
            return val
        if prop == "power_factor":
            val = bzana.get_power_factor(output=output)[dop_type][temp][index]
            return val
        if prop == "electronic_conductivity":
            val = bzana.get_conductivity(output=output)[dop_type][temp][index]
            return val
        if prop == "el_thermal_conductivity":
            val = bzana.get_thermal_conductivity(output=output)[dop_type][temp][index]
            return val
        if prop == "zt":
            val = bzana.get_zt(output=output)[dop_type][temp][index]
            return val
    else:
        print(
            "Value not available in the analyzer, run boltztrap with appropriate settings"
        )


if __name__ == "__main__":
    vrun = path = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-RELAX-bulk@mp_149/vasprun.xml",
        )
    )
    b = boltz_run(vrun)
    val = get_prop(b, prop="zt")
    print(val)
