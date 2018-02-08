# coding: utf-8
from math import pi
from pymatgen.core.structure import Structure
from operator import itemgetter
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
import os
from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator
import numpy as np
from math import log
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.core.periodic_table import Element
import json,sys
from pymatgen.core.structure import Structure
import numpy as np
import xml.etree.ElementTree as ET
from monty.serialization import loadfn

def get_descrp_arr(elm=''):
      arr=[]
      try:     
        tmp=str(os.path.join(os.path.dirname(__file__),'Elements.json'))
        f=open(tmp,'r')
        dat=json.load(f)
        f.close()
        d=dat[elm]
        arr=[]
        for k,v in d.items():
           arr.append(v)
        arr=np.array(arr).astype(float)
        #print os.path.join(os.path.dirname(__file__),'Elements.json')
      except:  
         pass
      return arr
def packing_fraction(s=None):
    total_rad = 0
    for site in s:
      total_rad += site.specie.atomic_radius ** 3
    pf=np.array([4 * pi * total_rad / (3 * s.volume)])
    return pf


def get_rdf(s=None,cutoff=10.0,intvl=0.1):
    neighbors_lst = s.get_all_neighbors(cutoff)
    all_distances = np.concatenate(tuple(map(lambda x: \
     [itemgetter(1)(e) for e in x], neighbors_lst)))
    rdf_dict = {}
    dist_hist, dist_bins = np.histogram(all_distances, \
     bins=np.arange(0, cutoff + intvl, intvl), density=False)
    shell_vol = 4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) \
    - np.power(dist_bins[:-1], 3))
    number_density = s.num_sites / s.volume
    rdf = dist_hist / shell_vol / number_density/len(neighbors_lst)
    return [{'distances': dist_bins[:-1], 'distribution': rdf}]


def get_comp_descp(struct=''):  
       cat=[]
       try: 
        if len(struct)<50: #TO AVOID SEG FAULT IN SYMMETRY LIBRARY E.G. mp-686203
          spg=str(SpacegroupAnalyzer(struct).get_space_group_number()).split()[0]
          s=(struct).get_primitive_structure()
        else:
          s=(struct)
          spg=1
        a=s.lattice.abc[0]
        b=s.lattice.abc[1]
        c=s.lattice.abc[2]
        ab=float(s.lattice.abc[0])/float(s.lattice.abc[1])
        bc=float(s.lattice.abc[1])/float(s.lattice.abc[2])
        ca=float(s.lattice.abc[2])/float(s.lattice.abc[0])
        alph=s.lattice.angles[0]
        bet=s.lattice.angles[1]
        gam=s.lattice.angles[2]
        #print 'spg=',spg
        v_pa=float(s.volume)/float(s.composition.num_atoms)
        logv_pa=log(float(s.volume)/float(s.composition.num_atoms))
        cell=np.array([a,b,c,ab,bc,ca,alph,bet,gam,v_pa,logv_pa,s.composition.num_atoms,spg])
        pf=packing_fraction(s)
        rdf=np.array(get_rdf(s=s)[0]['distribution'])
        #max_erdf=np.array(get_erdf(struct=s,val='max')[0]['distribution'])
        #min_erdf=np.array(get_erdf(struct=s,val='min')[0]['distribution'])
        
        comp=s.composition
        el_dict=comp.get_el_amt_dict()
        arr=[]
        #print ("eldi",el_dict,type(el_dict))
        for k,v in el_dict.items():
	    #print k,v
            des=get_descrp_arr(k)
	    #print k,v,des
            arr.append(des)
        #print arr,len(arr)
        #arr=np.array(arr)
        #print len(arr),arr
	#arr=np.concatenate((arr,cell),axis=0)
        #print 'arr=',arr
        mean_des=np.mean(arr,axis=0)
        max_des=np.maximum.reduce(arr)
        min_des=np.minimum.reduce(arr)
        sum_des=np.sum(arr,axis=0)
	#var_des=np.var(arr,axis=0)
	#diff_des=np.fabs(np.diff(arr,axis=0))[0]
        #print mean_des,len(mean_des)
	#cat=np.concatenate((mean_des,cell),axis=0).astype(float)
        cat=np.concatenate((mean_des,max_des,min_des,sum_des,cell,pf,rdf),axis=0).astype(float)

        #print cat,len(cat)

        #print 'org',cat,len(cat)
        ##cat=(cat)/float(np.mean(np.array(cat)))
        #print cat,len(cat)
	#cat=np.concatenate((mean_des,sum_des,std_des,cell),axis=0)
	#cat=np.concatenate((mean_des,sum_des,diff_des,std_des,cell),axis=0)
	#print mean_des
	#print sum_des
	#print diff_des,len(diff_des)
	#print std_des,len(std_des)
        #print cat,len(cat)
       except: 
        pass
       return cat
