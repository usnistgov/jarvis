"""
Classical Force-field Inspired Descriptors (CFID)
Find details in:
https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
"""
from __future__ import unicode_literals, print_function
from jarvis.core.atoms import Atoms
from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.core.specie import Specie
import matplotlib.pyplot as plt

plt.switch_backend("agg")
from collections import defaultdict
import itertools

# from scipy.stats import gaussian_kde
from math import pi
from operator import itemgetter
import collections, math, os
import numpy as np
from math import log
import json, sys
import numpy as np

# Note time limit for angular part is hardcoded

timelimit = 300  # 5 minutes


class CFID(object):
    def __init__(self, atoms):
        self._atoms = atoms

    def get_comp_descp(
        self,
        jcell=True,
        jmean_chem=True,
        jmean_chg=True,
        jrdf=False,
        jrdf_adf=True,
        print_names=False,
    ):
        """
        Get chemo-structural CFID decriptors

        Args:
        struct: Structure object
        jcell: whether to use cell-size descriptors
        jmean_chem: whether to use average chemical descriptors
        jmean_chg: whether to use average charge distribution descriptors
        jmean_rdf: whether to use radial distribution descriptors
        jrdf_adf: whether to use radial as well as angle distribution descriptors
        print_names: whether to print names of descriptors
        Returns:
          cat: catenated final descriptors
        """
        cat = []
        s = self._atoms
        cell = []
        mean_chem = []
        rdf = []
        adf = []
        nn = []
        mean_chg = []
        adfa = []
        adfb = []
        ddf = []

        if jmean_chem == True:
            el_dict = s.composition._content
            # print (el_dict,type(el_dict))
            arr = []
            for k, v in el_dict.items():
                des = Specie(k).get_descrp_arr
                arr.append(des)
            mean_chem = np.mean(arr, axis=0)
            # print ('mean_chem',len(mean_chem))

        if jcell == True:
            v_pa = round(float(s.volume) / float(s.num_atoms), 5)
            logv_pa = round(log(v_pa), 5)
            pf = s.packing_fraction
            density = round(s.density, 5)
            cell = np.array([v_pa, logv_pa, pf, density])
            # print ('jcell',len(cell))

        if jrdf == True:
            Nbrs = NeighborsAnalysis(s)
            _, distrdf, nn = Nbrs.get_rdf()
            rdf = np.array(distrdf)
            print("rdf", len(rdf))

        if jrdf_adf == True:
            try:
                adfa = np.zeros(179)
                adfb = np.zeros(179)
                ddf = np.zeros(179)
                rdf = np.zeros(100)
                nn = np.zeros(100)
                distributions = NeighborsAnalysis(s).get_all_distributions
                rdf = distributions["rdf"]
                nn = distributions["nn"]
                adfa = distributions["adfa"]
                adfb = distributions["adfb"]
                ddf = distributions["ddf"]

            except:
                pass
            adfa = np.array(adfa)
            adfb = np.array(adfb)
            rdf = np.array(rdf)
            ddf = np.array(ddf)
            nn = np.array(nn)
            # print ('adfa',len(adfa))
            # print ('ddf',len(ddf))
            # print ('adfb',len(adfb))
            # print ('rdf',len(rdf))
            # print ('nn',len(nn))

        if jmean_chg == True:
            chgarr = []
            el_dict = s.composition._content
            for k, v in el_dict.items():
                chg = Specie(k).get_chgdescrp_arr
                chgarr.append(chg)
            mean_chg = np.mean(chgarr, axis=0)
            # print ('mean_chg',len(mean_chg))

        if print_names == True:
            nmes = []
            chem_nm = get_descrp_arr_name()
            for d, nm in zip(
                [mean_chem, cell, mean_chg, rdf, adfa, adfb, ddf, nn],
                ["mean_chem", "cell", "mean_chg", "rdf", "adfa", "adfb", "ddf", "nn"],
            ):
                if d != []:
                    for ff, dd in enumerate(d):
                        cat.append(dd)
                        if nm == "mean_chem":
                            tag = chem_nm[ff]
                        else:
                            tag = str(nm) + str("_") + str(ff)
                        nmes.append(str(tag))
            cat = np.array(cat).astype(float)
            # print (nmes,len(nmes))
            return nmes
        else:
            for d, nm in zip(
                [mean_chem, cell, mean_chg, rdf, adfa, adfb, ddf, nn],
                ["mean_chem", "cell", "mean_chg", "rdf", "adfa", "adfb", "ddf", "nn"],
            ):
                if d != []:
                    for ff, dd in enumerate(d):
                        cat.append(dd)
            cat = np.array(cat).astype(float)
        return cat


if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    cfid = CFID(Si).get_comp_descp()
    print(len(cfid))
