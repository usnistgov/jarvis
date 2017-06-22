from sklearn.neural_network import MLPClassifier
from pymatgen.core.lattice import Lattice
from ase import data
import numpy as np
from asap3.analysis.rdf import RadialDistributionFunction
from pymatgen.io.aseio import AseAtomsAdaptor
from pymatgen.core.structure import Structure
import json,pickle

def YuriStructure(strt=None):
    """
    Takes Pymatgen structure and gives structure in Yuri's format
    """
    info={}
    latt=np.array([strt.lattice.matrix])
    info['lattice']=latt.tolist()
    arr=[]
    for site in strt:
        arr.append([str(site.specie),site.coords[0],site.coords[1],site.coords[2]])
    info['sites']=arr
    return info

def pmg_yuri(yuri):
    sites=yuri['sites']
    natoms=len(sites)
    x=np.zeros((natoms))
    y=np.zeros((natoms))
    z=np.zeros((natoms))
    lat=Lattice(yuri['lattice'])
    typ= np.empty((natoms),dtype="S20")
    coords = list()#np.zeros((natoms))
    for i, j in enumerate(sites):
                 typ[i]=(j[0])
                 x[i]=j[1]
                 y[i]=j[2]
                 z[i]=j[3]
                 coords.append([x[i],y[i],z[i]])
    struct=Structure(lat,typ,coords,coords_are_cartesian=True)
    return struct


def expand_c(strt=None,encut=500,length=50,minus=-10.0,plus=10.0,intvl=20):
    arr=[]
    for i,percent in enumerate(np.linspace(minus,plus,intvl)):
        strain=float(percent)/float(100.0)
        tmp=strt.copy()
        tmp.apply_strain(strain)
        #pos=Poscar(tmp)
        #pos.comment=str('Percentage-')+str(percent)
        arr.append(tmp)
    return arr
def rdf(contcar):
        rng=15.0
        bins = 200
        info={}
        atoms= AseAtomsAdaptor().get_atoms(contcar)
        try:
           sa1=int(float(rng)/float(max(abs(atoms.get_cell()[0]))))+2
           sa2=int(float(rng)/float(max(abs(atoms.get_cell()[1]))))+2
           sa3=int(float(rng)/float(max(abs(atoms.get_cell()[2]))))+2
           #print ('scal=',sa1,sa2,sa3)
           #atoms=atoms*(int(float(rng)/float(atoms.get_cell()[0][0]))+2,int(float(rng)/float(atoms.get_cell()[1][1]))+2,int(float(rng)/float(atoms.get_cell()[2][2]))+2)
           atoms=atoms*(sa1,sa2,sa3)
        except:
           atoms=atoms*(9,9,9)
        symb=atoms.get_chemical_symbols()
        symbs=[]
        for a in symb:
            if a not in symbs:
               symbs.append(a)
        symb_comb=[]
        x = np.arange(bins) * rng / bins
        RDFobj = RadialDistributionFunction(atoms, rng, bins, verbose=True)
        symb_comb=[]
        #print ("symbs",symbs)
        #plt.close()
        if len(symbs)>1:
           for el1 in  symbs:
              for el2 in  symbs:
                 el=str(el1)+str('-')+str(el2)
                 if str(el1)+str('-')+str(el2) not in symb_comb and str(el2)+str('-')+str(el1) not in symb_comb:
                     symb_comb.append(el)
                     rdf = RDFobj.get_rdf(elements=(data.atomic_numbers[str(el1)], data.atomic_numbers[str(el2)]))
                     #plt.tight_layout()
                     #plt.plot(x, rdf,label=str(el1)+str('-')+str(el2),linewidth=2)
                     key=str(el1)+str('-')+str(el2)
                     value=rdf
                     info[key]=value
        elif len(symbs)==1:
                     rdf = RDFobj.get_rdf(elements=(data.atomic_numbers[str(symbs[0])], data.atomic_numbers[str(symbs[0])]))
                     key=str(symbs[0])+str('-')+str(symbs[0])
                     value=rdf
                     info[key]=value
        return info


#s=Structure.from_file('POSCAR')
#ec=expand_c(s)

#for i in ec:

#   l=rdf(i)
#   print l

def train(filename='NN3_data.json'):
  f=open(filename,'r')
  dat=json.load(f)
  f.close()
  inv_arr=[]
  eng_arr=[]
  for i in dat:
    for j in i['info']:
        k=pmg_yuri(j['struct'])
        l=rdf(k)
        e=j['energy']
        nat=len(j['struct']['sites'])
        e_at=float(e)/float(nat)
        #print l['Si-Si']
        inv_arr.append(l['Si-Si'])
        eng_arr.append(e_at)
        #print e_at
  eng_arr = np.asarray(eng_arr, dtype="|S32")
  clf = MLPClassifier(solver='lbfgs', alpha=1e-14,hidden_layer_sizes=(37, 37), random_state=1)
  clf.fit(inv_arr,eng_arr)
  pickle.dump( clf, open( "save.p", "wb" ) )


def pred(file1='NN3_data.json',file2="save.p"):
  object_file =pickle.load( open( file2, "rb" ) )
  f=open(file1,'r')
  dat=json.load(f)
  f.close()
  inv_arr=[]
  eng_arr=[]
  for i in dat:
    for j in i['info']:
        k=pmg_yuri(j['struct'])
        l=rdf(k)
        e=j['energy']
        nat=len(j['struct']['sites'])
        e_at=float(e)/float(nat)
        inv=l['Si-Si']
        print object_file.predict(inv),e_at
pred()
