from jarvis.core.atoms import Atoms
import numpy as np
import xmltodict
from collections import OrderedDict
from jarvis.core.atoms import Atoms
import numpy as np
import xmltodict
from collections import OrderedDict
from jarvis.core.kpoints import Kpoints3D as Kpoints


class TotalDos(object):
    def __init__(
        self,
        energies=[],
        values=[],
        integrated_values=[],
        spin=None,
        element=None,
        orbital=None,
    ):
        self.energues = energies
        self.values = values
        self.integrated_values = integrated_values
        self.element = element
        self.orbital = orbital
        self.spin = spin


class Vasprun(object):
    def __init__(self, filename="vasprun.xml", data={}):
        self._filename = filename
        self._data = data
        self.ionic_steps = None
        self.electronic_steps = None
        self.input_parameters = None
        if self._data == {}:
            self.xml_to_dict()

    def xml_to_dict(self):
        with open(self._filename) as fd:
            data = xmltodict.parse(fd.read())
            self._data = data
            self.ionic_steps = data["modeling"]["calculation"]
            if type(self.ionic_steps) is not list:
                self.ionic_steps = [self.ionic_steps]
            if self.input_parameters is None:
                self.input_parameters = self.all_input_parameters

    @property
    def final_energy(self):
        return float(self.ionic_steps[-1]["scstep"][-1]["energy"]["i"][11]["#text"])

    @property
    def efermi(self):
        return float(self.ionic_steps[-1]["dos"]["i"]["#text"])

    @property
    def num_atoms(self):
        return self._data["modeling"]["atominfo"]["atoms"]

    @property
    def num_types(self):
        return int(self._data["modeling"]["atominfo"]["types"])

    @property
    def elements(self):
        elements = [
            self._data["modeling"]["atominfo"]["array"][0]["set"]["rc"][i]["c"][0]
            for i in range(
                len(self._data["modeling"]["atominfo"]["array"][0]["set"]["rc"])
            )
        ]
        if len(elements) != self.num_atoms:
            ValueError("Number of atoms is  not equal to number of elements")
        elements = [str(i) for i in elements]
        return elements

    def vrun_structure_to_atoms(self, s={}):
        lattice_mat = np.array(
            [[float(j) for j in i.split()] for i in s["crystal"]["varray"][0]["v"]]
        )
        frac_coords = np.array(
            [[float(j) for j in i.split()] for i in s["varray"]["v"]]
        )
        elements = self.elements
        atoms = Atoms(
            lattice_mat=lattice_mat,
            elements=elements,
            coords=frac_coords,
            cartesian=False,
        )
        return atoms

    @property
    def all_energies(self):
        energies = []
        for i in self.ionic_steps:
            en = float(i["energy"]["i"][1]["#text"])
            energies.append(en)
        return np.array(energies)

    @property
    def is_spin_polarized(self):
        if self.all_input_parameters["ISPIN"] == "2":
            return True
        else:
            return False

    @property
    def all_structures(self):
        structs = []
        for i in self.ionic_steps:
            s = i["structure"]
            atoms = self.vrun_structure_to_atoms(s)
            structs.append(atoms)
        return structs

    @property
    def eigenvalues(self):
        nkpts = len(self.kpoints._kpoints)
        all_up_eigs = []
        all_dn_eigs = []

        for j in range(nkpts):
            eigs = np.array(
                [
                    [float(jj) for jj in ii.split()]
                    for ii in (
                        self.ionic_steps[-1]["eigenvalues"]["array"]["set"]["set"][0]
                    )["set"][j]["r"]
                ]
            )
            all_up_eigs.append(eigs)
        for j in range(nkpts):
            eigs = np.array(
                [
                    [float(jj) for jj in ii.split()]
                    for ii in (
                        self.ionic_steps[-1]["eigenvalues"]["array"]["set"]["set"][1]
                    )["set"][j]["r"]
                ]
            )
            all_dn_eigs.append(eigs)
        all_up_eigs = np.array(all_up_eigs)
        all_dn_eigs = np.array(all_dn_eigs)
        return all_up_eigs, all_dn_eigs

    @property
    def all_forces(self):
        forces = []
        for m in self.ionic_steps:
            force = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][0]["v"]]
            )

            forces.append(force)
        return np.array(forces)

    @property
    def all_stresses(self):
        stresses = []
        for m in self.ionic_steps:
            stress = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][1]["v"]]
            )

            stresses.append(stress)
        return np.array(stresses)

    @property
    def all_input_parameters(self):
        d = OrderedDict()
        # import type
        for i in self._data["modeling"]["parameters"]["separator"]:
            for j, k in i.items():
                if j == "i":
                    for m in k:
                        if "#text" in m:
                            d[m["@name"]] = m["#text"]
                else:
                    if type(k) is list:
                        for n in k:
                            for p, q in n.items():
                                if p == "i":
                                    for r in q:
                                        if "#text" in r:
                                            d[r["@name"]] = r["#text"]
                                else:
                                    if type(q) is list:
                                        for s in q:
                                            if "#text" in s:
                                                d[s["@name"]] = s["#text"]
        return d

    @property
    def kpoints(self):
        kplist = np.array(
            [
                [float(j) for j in i.split()]
                for i in self._data["modeling"]["kpoints"]["varray"][0]["v"]
            ]
        )
        kpwt = np.array(
            [float(i) for i in self._data["modeling"]["kpoints"]["varray"][1]["v"]]
        )
        return Kpoints(kpoints=kplist, kpoints_weights=kpwt)

    def get_bandstructure(self, spin=0):
        plt.clf()
        for i, ii in enumerate(self.eigenvalues[spin][:, :, 0].T - float(v.efermi)):
            plt.plot(ii, color="b")
        return plt

    @property
    def total_dos(self):
        energies = []
        spin_up = []
        spin_up_data = np.array(
            [
                [float(j) for j in i.split()]
                for i in self.ionic_steps[-1]["dos"]["total"]["array"]["set"]["set"][0][
                    "r"
                ]
            ]
        )
        energies = spin_up_data[:, 0]
        spin_up = spin_up_data[:, 1]
        if self.is_spin_polarized:
            spin_dn = []
            spin_dn_data = np.array(
                [
                    [float(j) for j in i.split()]
                    for i in self.ionic_steps[-1]["dos"]["total"]["array"]["set"][
                        "set"
                    ][1]["r"]
                ]
            )
            spin_dn = -1 * spin_dn_data[:, 1]
        return energies, spin_up, spin_dn


class Oszicar(object):
     def __init__(self,filename, data={}):
          self.filename= filename
          self.data = data
          if self.data=={}:
             f=open(filename,'r')
             lines=f.read().splitlines()
             f.close()
             self.data = lines

     def magnetic_moment(self):
         return self.ionic_steps()[-1][-1]

     def ionic_steps(self):
         ionic_data = []
         for i in self.data:
           if 'E0' in i:
             ionic_data.append(i.split())
         return ionic_data
               
         
class Outcar(object):
      def __init__(self,filename,data={}):
          self.filename= filename
          self.data = data
          if self.data=={}:
             f=open(filename,'r')
             lines=f.read().splitlines()
             f.close()
             self.data = lines
        
      
class Chgcar(object):
      def __init__(self,filename,data={}):
          self.filename= filename
          self.data = data
          if self.data=={}:
             f=open(filename,'r')
             lines=f.read().splitlines()
             f.close()
             self.data = lines


if __name__ == "__main__":
    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/OSZICAR"
    osz=Oszicar(filename=filename)
    print (osz.ionic_steps(),osz.magnetic_moment())
    import sys
    sys.exit()


    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/vasprun.xml"
    v = Vasprun(filename=filename)
    # print (v._filename, v.final_energy)
    # print (v._filename,'elements', v.elements)
    print(v._filename, "structs", v.all_structures)
    # print ('pmg',Vasprun(v._filename).final_energy)
    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-BAND-bulk@mp_541837/vasprun.xml"
    v = Vasprun(filename=filename)
    print(v._filename, "structs", v.all_structures)
    # print (v._filename,v.final_energy)
    # print (v._filename,v.elements)
    # print ('pmg',Vasprun(v._filename).final_energy)
    # print (v._data.keys())
    # print ()
    # print ()
    # print (v._data['modeling'].keys())
    # print ()
    # print ()
    ##'generator', 'incar', 'kpoints', 'parameters', 'atominfo', 'structure', 'calculation'
    # print ('incar',v._data['modeling']['incar'])
    # print ()
    # print ()
    # print ('kpoints',v._data['modeling']['kpoints'])
    # print ()
    # print ()
    # print ('atominfo',v._data['modeling']['atominfo'])
    # print ()
    # print ()
    # print ('parameters',v._data['modeling']['parameters'])
    # print ()
    # print ()
    # print ('structure',v._data['modeling']['structure'])
    # print ()
    # print ()
    # print ('calculation',v._data['modeling']['calculation'].keys())
    # print ()
    # print ()
    ##calculation odict_keys(['scstep', 'structure', 'varray', 'energy', 'time', 'eigenvalues', 'separator', 'dos', 'projected'])
