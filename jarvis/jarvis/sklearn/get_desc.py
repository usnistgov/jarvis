"""
Classical Force-field Inspired Descriptors (CFID)
Find details in:
https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
"""
from __future__ import  unicode_literals, print_function
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from collections import defaultdict
import itertools
from scipy.stats import gaussian_kde
from math import pi
from operator import itemgetter
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
#from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
import collections,math,os
#from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator
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
from pandas import DataFrame as pd
from pymatgen.core.lattice import Lattice
from pymatgen.core import Composition

el_chrg_json=str(os.path.join(os.path.dirname(__file__),'element_charge.json')) 
el_chem_json=str(os.path.join(os.path.dirname(__file__),'Elements.json')) 

def get_effective_structure(s=None,tol=8.0):
  
  coords=s.cart_coords
  range_x=max(coords[:,0])-min(coords[:,0])
  range_y=max(coords[:,1])-min(coords[:,1])
  range_z=max(coords[:,2])-min(coords[:,2])
  a=s.lattice.matrix[0][0]
  b=s.lattice.matrix[1][1]
  c=s.lattice.matrix[2][2]
  #print 'range_x=max',range_x
  if abs(a-range_x)>tol:a=range_x+tol
  if abs(b-range_y)>tol:b=range_y+tol
  if abs(c-range_z)>tol:c=range_z+tol
  #print 'abc',a,b,c
  arr=Lattice([[a,s.lattice.matrix[0][1],s.lattice.matrix[0][2]],[s.lattice.matrix[1][0],b,s.lattice.matrix[1][2]],[s.lattice.matrix[2][0],s.lattice.matrix[2][1],c]])
  s=Structure(arr,s.species,coords,coords_are_cartesian=True)
  return s

def el_combs(s=[]):
  symb=s.symbol_set
  tmp=map('-'.join, itertools.product(symb, repeat=2))
  comb=list(set([str('-'.join(sorted(i.split('-')))) for i in tmp]))
  return comb

def flatten_out(arr=[],tol=0.1):
    rcut_buffer=tol 
    
    io1=0
    io2=1
    io3=2
    delta=arr[io2]-arr[io1]
     
    while delta < rcut_buffer and io3 <len(arr):
      io1=io1+1
      io2=io2+1
      io3=io3+1
      delta=arr[io2]-arr[io1]
     # print ('arr',len(arr),io1,io2,io3)
    rcut=(arr[io2]+arr[io1])/float(2.0)
    return rcut

def smooth_kde(x,y):
    denn=gaussian_kde(y)
    denn.covariance_factor = lambda : .25
    denn._compute_covariance()
    xs = np.linspace(0,max(x),100)
    kde=denn(xs)
    return kde

def get_prdf(s=None,comb='',cutoff=10.0,intvl=0.1,rdf_tol=0.1,plot_prdf=False,filename='prdf.png'):
    
    neighbors_lst = s.get_all_neighbors(cutoff)
    comb=el_combs(s=s)   #[str(i[0])+str('-')+str(i[1]) for i in list(itertools.product(sps, repeat=2))]
    info={}
    for c in comb:
     for i,ii in enumerate(neighbors_lst):
      for j in ii:
         comb1=str(s[i].specie)+str('-')+str( j[0].specie)
         comb2=str(j[0].specie)+str('-')+str( s[i].specie)
         if comb1==c or comb2==c:
           info.setdefault(c, []).append(j[1])
    plt.rcParams.update({'font.size': 22})
    for i in info.items():
      i[1].sort()
      dist_hist, dist_bins = np.histogram(i[1], bins=np.arange(0, cutoff + intvl, intvl), density=False)
      shell_vol = 4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
      number_density = s.num_sites / s.volume
      rdf = dist_hist / shell_vol / number_density/len(neighbors_lst)
      #newy=smooth_kde(dist_bins[:-1],rdf)
      if plot_prdf==True:
      #plt.plot(newy,rdf,label=i[0],linewidth=2)
        plt.plot(dist_bins[:-1],rdf,label=i[0],linewidth=2)
    
    if plot_prdf==True:

     plt.legend(prop={'size':26})
     plt.xlabel('r ($\AA$)')
     plt.ylabel('g(r)')
     plt.ylim(ymin=0)
     plt.tight_layout()
     plt.savefig(filename)
     plt.close()

    cut_off={}

    for i,j in info.items():
        cut_off[i]=flatten_out(arr=j,tol=0.1)
   # print 'cut_off',cut_off
   # for i,j in info.items():
   #   print
   #   print i,sorted(set(j))
   #   print
    return max(cut_off.items(), key=itemgetter(1))[1]


def get_rdf(s=None,cutoff=10.0,intvl=0.1):
    neighbors_lst = s.get_all_neighbors(cutoff)
    all_distances = np.concatenate(tuple(map(lambda x: \
     [itemgetter(1)(e) for e in x], neighbors_lst)))
    rdf_dict = {}
    dist_hist, dist_bins = np.histogram(all_distances, \
     bins=np.arange(0, cutoff + intvl, intvl), density=False) #equivalent to bon-order
    shell_vol = 4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) \
    - np.power(dist_bins[:-1], 3))
    number_density = s.num_sites / s.volume
    rdf = dist_hist / shell_vol / number_density/len(neighbors_lst)
    return dist_bins[:-1],[round(i,4) for i in rdf],dist_hist/float(len(s)) #[{'distances': dist_bins[:-1], 'distribution': rdf}]
    #bins,rdf,nearest neighbour

def rdf_ang_dist(s='',c_size=10.0,plot=True,max_cut=5.0):
    x,y,z=get_rdf(s)
    arr=[]
    for i,j in zip(x,z):
      if j>0.0:
        arr.append(i)
    box=s.lattice.matrix
    rcut_buffer=0.11
    io1=0
    io2=1
    io3=2
    #print 'here1'
    delta=arr[io2]-arr[io1]
    #while (delta < rcut_buffer and io2<len(arr)-2):
    while (delta < rcut_buffer and arr[io2]<max_cut):
      io1=io1+1
      io2=io2+1
      io3=io3+1
      delta=arr[io2]-arr[io1]
    #print 'here2'
    rcut1=(arr[io2]+arr[io1])/float(2.0)
    #print "arr=",arr[0],arr[1],arr[2],arr[3]
    rcut=get_prdf(s=s) # (arr[io2]+arr[io1])/float(2.0)
    #print "rcut=",rcut,"rcut1=",rcut1
    #print io3,io2,len(arr)
    delta=arr[io3]-arr[io2]
    while (delta < rcut_buffer and arr[io3]<max_cut and arr[io2]<max_cut):
      io2=io2+1
      io3=io3+1
      #print arr[io3],arr[io2],io3,io2
      #print arr[io3],arr[io2]
      delta=arr[io3]-arr[io2]
    rcut2=float(arr[io3]+arr[io2])/float(2.0)
    #print "rcut2=",rcut2

    #print 'here3'
    dim1=int(float(c_size)/float( max(abs(box[0]))))+1
    dim2=int(float(c_size)/float( max(abs(box[1]))))+1
    dim3=int(float(c_size)/float( max(abs(box[2]))))+1
    dim=[dim1,dim2,dim3]
    #rcut2=(arr[2]+arr[3])/2.0
    dim=np.array(dim)
    coords= s.frac_coords #get_frac(s.cart_coords,s.lattice.matrix) #s.frac_coords

    lat=np.zeros((3,3))
    lat[0][0]=dim[0]*box[0][0]
    lat[0][1]=dim[0]*box[0][1]
    lat[0][2]=dim[0]*box[0][2]
    lat[1][0]=dim[1]*box[1][0]
    lat[1][1]=dim[1]*box[1][1]
    lat[1][2]=dim[1]*box[1][2]
    lat[2][0]=dim[2]*box[2][0]
    lat[2][1]=dim[2]*box[2][1]
    lat[2][2]=dim[2]*box[2][2]
    #print "lat="
    #print lat
    #print ""  

    all_symbs=[i.symbol for i in s.species]
    nat=len(coords)
    new_nat=nat*dim[0]*dim[1]*dim[2]
    new_coords=np.zeros((new_nat,3))
    new_symbs= [] #np.chararray((new_nat))

    count=0
    for i in range(nat):
      for j in range(dim[0]):
       for k in range(dim[1]):
        for l in range(dim[2]):
          new_coords[count][0]=(coords[i][0]+j)/float(dim[0])
          new_coords[count][1]=(coords[i][1]+k)/float(dim[1])
          new_coords[count][2]=(coords[i][2]+l)/float(dim[2])
          new_symbs.append(all_symbs[i])
          count=count+1

    #print 'here4'
    nat=new_nat
    coords=new_coords
    znm=0
    bond_arr=[]
    deg_arr=[]
    nn=np.zeros((nat),dtype='int')
    max_n=500 #maximum number of neighbors
    dist=np.zeros((max_n,nat))
    nn_id=np.zeros((max_n,nat),dtype='int')
    bondx=np.zeros((max_n,nat))
    bondy=np.zeros((max_n,nat))
    bondz=np.zeros((max_n,nat))
    dim05=[float(1/2.) for i in dim]
    #print "dim",dim
    #print "dim05",dim05
    for i in range(nat):
     for j in range(i+1,nat):
       diff=coords[i]-coords[j]
       for v in range(3):
         if np.fabs(diff[v])>=dim05[v]:
           diff[v]=diff[v]-np.sign(diff[v])
       new_diff=np.dot(diff,lat)
       dd=np.linalg.norm(new_diff)
       if dd<rcut and dd>=0.1:
           nn_index=nn[i] #index of the neighbor
           nn[i]=nn[i]+1
           dist[nn_index][i]=dd #nn_index counter id
           nn_id[nn_index][i]=j #exact id
           bondx[nn_index][i]=new_diff[0]
           bondy[nn_index][i]=new_diff[1]
           bondz[nn_index][i]=new_diff[2]
           nn_index1=nn[j] #index of the neighbor
           nn[j]=nn[j]+1
           dist[nn_index1][j]=dd #nn_index counter id
           nn_id[nn_index1][j]=i #exact id
           bondx[nn_index1][j]=-new_diff[0]
           bondy[nn_index1][j]=-new_diff[1]
           bondz[nn_index1][j]=-new_diff[2]


    ang_at={}

    for i in range(nat):
     for in1 in range(nn[i]):
       j1=nn_id[in1][i]
       for in2 in range(in1+1,nn[i]):
         j2=nn_id[in2][i]
         nm=dist[in1][i]*dist[in2][i]
         if nm!=0:
          rrx=bondx[in1][i]*bondx[in2][i]
          rry=bondy[in1][i]*bondy[in2][i]
          rrz=bondz[in1][i]*bondz[in2][i]
          cos=float(rrx+rry+rrz)/float(nm)
          if cos<=-1.0:
           cos=cos+0.000001
          if cos>=1.0:
           cos=cos-0.000001
          deg=math.degrees(math.acos(cos))
          #ang_at.setdefault(deg, []).append(i)
          ang_at.setdefault(round(deg,3), []).append(i)
         else:
           znm=znm+1
    angs=np.array([float(i) for i in ang_at.keys()])
    #print "Angs",angs
    norm=np.array([float(len(i))/float(len(set(i))) for i in ang_at.values()])
    #print ('angs',angs,len(angs))
    #angs_np=np.array(angs)
    #print ('angs',angs_np,len(angs_np))
    #print ('norm',norm.shape,angs.shape) 
    ang_hist1, ang_bins1 = np.histogram(angs,weights=norm,bins=np.arange(1, 181.0, 1), density=False)

    #print 'here5'

    #print "rcut1=",rcut1
    # Dihedral angle distribution
    znm=0
    bond_arr=[]
    deg_arr=[]
    nn=np.zeros((nat),dtype='int')
    max_n=500 #maximum number of neighbors
    dist=np.zeros((max_n,nat))
    nn_id=np.zeros((max_n,nat),dtype='int')
    bondx=np.zeros((max_n,nat))
    bondy=np.zeros((max_n,nat))
    bondz=np.zeros((max_n,nat))
    dim05=[float(1/2.) for i in dim]
    for i in range(nat):
     for j in range(i+1,nat):
       diff=coords[i]-coords[j]
       for v in range(3):
         if np.fabs(diff[v])>=dim05[v]:
           diff[v]=diff[v]-np.sign(diff[v])
       new_diff=np.dot(diff,lat)
       dd=np.linalg.norm(new_diff)
       if dd<rcut1 and dd>=0.1:
           nn_index=nn[i] #index of the neighbor
           nn[i]=nn[i]+1
           #print nn_index
           dist[nn_index][i]=dd #nn_index counter id
           nn_id[nn_index][i]=j #exact id
           bondx[nn_index,i]=new_diff[0]
           bondy[nn_index,i]=new_diff[1]
           bondz[nn_index,i]=new_diff[2]
           nn_index1=nn[j] #index of the neighbor
           nn[j]=nn[j]+1
           dist[nn_index1][j]=dd #nn_index counter id
           nn_id[nn_index1][j]=i #exact id
           bondx[nn_index1,j]=-new_diff[0]
           bondy[nn_index1,j]=-new_diff[1]
           bondz[nn_index1,j]=-new_diff[2]
    #print 'here6',rcut2

    dih_at={}
    for i in range(nat):
     for in1 in range(nn[i]):
     #for in1 in range(1):
       j1=nn_id[in1][i]
       if (j1 > i):
        """
        # angles between i,j, k=nn(i), l=nn(i)
        for in2 in range(nn[i]):    # all other nn of i that are not j
            j2=nn_id[in2][i]
            if (j2 != j1):
               for in3 in range(in2+1,nn[i]):  
                   j3=nn_id[in3][i]
                   if (j3 != j1):
                      
        # angles between i,j, k=nn(j), l=nn(j)
        for in2 in range(nn[j1]):    # all other nn of j that are not i
            j2=nn_id[in2][j1]
            if (j2 != i):
               for in3 in range(in2+1,nn[j1]):    
                   j3=nn_id[in3][j1]
                   if (j3 != i):
        """
        # angles between i,j, k=nn(i), l=nn(j)
        for in2 in range(nn[i]):    # all other nn of i that are not j
            j2=nn_id[in2][i]
            if (j2 != j1):          
               for in3 in range(nn[j1]):    # all other nn of j that are not i
                   j3=nn_id[in3][j1]
                   if (j3 != i):
                      v1=[] 
                      v2=[]  
                      v3=[]  
                      v1.append(bondx[in2][i])
                      v1.append(bondy[in2][i])
                      v1.append(bondz[in2][i])
                      v2.append(-bondx[in1][i])
                      v2.append(-bondy[in1][i])
                      v2.append(-bondz[in1][i])
                      v3.append(-bondx[in3][j1])
                      v3.append(-bondy[in3][j1])
                      v3.append(-bondz[in3][j1])
                      v23 = np.cross(v2, v3)
                      v12 = np.cross(v1, v2)
                      theta = math.degrees(math.atan2(np.linalg.norm(v2)*np.dot(v1, v23),np.dot(v12, v23)))
                      if theta < 0.00001: theta = - theta
                      #print "theta=",theta
                      dih_at.setdefault(round(theta,3), []).append(i) 
    dih=np.array([float(i) for i in dih_at.keys()])
    dih1=set(dih)
    #print "dih",dih1
    norm=np.array([float(len(i))/float(len(set(i))) for i in dih_at.values()])
    
    dih_hist1, dih_bins1 = np.histogram(dih,weights=norm,bins=np.arange(1, 181.0, 1), density=False)
    

### 2nd neighbors 
    znm=0
    bond_arr=[]
    deg_arr=[]
    nn=np.zeros((nat),dtype='int')
    max_n=250 #maximum number of neighbors
    dist=np.zeros((max_n,nat))
    nn_id=np.zeros((max_n,nat),dtype='int')
    bondx=np.zeros((max_n,nat))
    bondy=np.zeros((max_n,nat))
    bondz=np.zeros((max_n,nat))
    dim05=[float(1/2.) for i in dim]
    for i in range(nat):
     for j in range(i+1,nat):
       diff=coords[i]-coords[j]
       for v in range(3):
         if np.fabs(diff[v])>=dim05[v]:
           diff[v]=diff[v]-np.sign(diff[v])
       new_diff=np.dot(diff,lat)
       dd=np.linalg.norm(new_diff)
       if dd<rcut2 and dd>=0.1:
           nn_index=nn[i] #index of the neighbor
           nn[i]=nn[i]+1
           dist[nn_index][i]=dd #nn_index counter id
           nn_id[nn_index][i]=j #exact id
           bondx[nn_index,i]=new_diff[0]
           bondy[nn_index,i]=new_diff[1]
           bondz[nn_index,i]=new_diff[2]
           nn_index1=nn[j] #index of the neighbor
           nn[j]=nn[j]+1
           dist[nn_index1][j]=dd #nn_index counter id
           nn_id[nn_index1][j]=i #exact id
           bondx[nn_index1,j]=-new_diff[0]
           bondy[nn_index1,j]=-new_diff[1]
           bondz[nn_index1,j]=-new_diff[2]


    ang_at={}

    for i in range(nat):
     for in1 in range(nn[i]):
       j1=nn_id[in1][i]
       for in2 in range(in1+1,nn[i]):
         j2=nn_id[in2][i]
         nm=dist[in1][i]*dist[in2][i]
         if nm!=0:
          rrx=bondx[in1][i]*bondx[in2][i]
          rry=bondy[in1][i]*bondy[in2][i]
          rrz=bondz[in1][i]*bondz[in2][i]
          cos=float(rrx+rry+rrz)/float(nm)
          if cos<=-1.0:
           cos=cos+0.000001
          if cos>=1.0:
           cos=cos-0.000001
          deg=math.degrees(math.acos(cos))
          #ang_at.setdefault(deg, []).append(i)
          ang_at.setdefault(round(deg,3), []).append(i)
         else:
           znm=znm+1
    angs=np.array([float(i) for i in ang_at.keys()])
    norm=np.array([float(len(i))/float(len(set(i))) for i in ang_at.values()])
 

    ang_hist2, ang_bins2 = np.histogram(angs,weights=norm,bins=np.arange(1, 181.0, 1), density=False)

    if plot==True:
     plt.plot(ang_bins2[:-1],ang_hist2)
     plt.savefig('ang.png')
     plt.close()
     plt.plot(x,y)
     plt.savefig('angrdf.png')
     plt.close()

    return ang_hist1,ang_hist2,dih_hist1,y,z #adfa,adfb,ddf,rdf,bondo


def get_chgdescrp_arr(elm=''):
      arr=[]
      try:      
        f=open(el_chrg_json,'r')
        emdat=json.load(f)
        f.close()
        arr=emdat[elm][0][1]
      except:  
         pass
      return arr


def get_descrp_arr_name(elm='Al'):
      arr=[]
      try:      
        f=open(el_chem_json,'r')
        dat=json.load(f)
        f.close()
       
        d=dat[elm]
        arr=[]
        for k,v in d.items():
            arr.append(k)
      except:  
          pass 
      return arr

def get_descrp_arr(elm=''):
      arr=[]
      try:      
        f=open(el_chem_json,'r')
        dat=json.load(f)
        f.close()
       
        d=dat[elm]
        arr=[]
        for k,v in d.items():
           arr.append(v)
        arr=np.array(arr).astype(float)
      except:  
        pass
      return arr

def packing_fraction(s=None):
    total_rad = 0
    for site in s:
       total_rad += site.specie.atomic_radius ** 3
    pf=np.array([4 * pi * total_rad / (3 * s.volume)])
    return pf[0]



def get_comp_descp(struct='',jcell=True,jmean_chem=True,jmean_chg=True,jrdf=False,jrdf_adf=True,print_names=False):  
        cat=[]
     # try: 
        s= get_effective_structure(struct)
       #	 print (len(s))
        cell=[]
        mean_chem=[]
        rdf=[]
        adf=[]
        nn=[]
        mean_chg=[]
        adfa=[]
        adfb=[]
        ddf=[]


        if jmean_chem==True:
         comp=s.composition
         el_dict=comp.get_el_amt_dict()
         arr=[]
         for k,v in el_dict.items():
            des=get_descrp_arr(k)
            arr.append(des)
         mean_chem=np.mean(arr,axis=0)
         #print ('mean_chem',len(mean_chem))



        if jcell==True:
         v_pa=round(float(s.volume)/float(s.composition.num_atoms),5)
         logv_pa=round(log(float(s.volume)/float(s.composition.num_atoms)),5)
         pf=round(packing_fraction(s),5)
         density=round(s.density,5)
         cell=np.array([v_pa,logv_pa,pf,density])
         #print ('jcell',len(cell))


        if jrdf==True:
         distrdf,bins,bo=get_rdf(s=s)
         rdf=np.array(distrdf)
         #print ('rdf',len(rdf))

        if jrdf_adf==True:
         adfa,adfb,ddf,rdf,nn=rdf_ang_dist(s=s,plot=False)
         adfa=np.array(adfa)
         adfb=np.array(adfb)
         rdf=np.array(rdf)
         ddf=np.array(ddf)
         nn=np.array(nn)
         #print ('adfa',len(adfa))
         #print ('ddf',len(ddf))
         #print ('adfb',len(adfb))
         #print ('rdf',len(rdf))
         #print ('nn',len(nn))



        if jmean_chg==True:
  
         chgarr=[]
         for k,v in el_dict.items():
            chg=get_chgdescrp_arr(k)
            chgarr.append(chg)
         mean_chg=np.mean(chgarr,axis=0)
         #print ('mean_chg',len(mean_chg))


        if print_names==True:
         nmes=[]
         chem_nm=get_descrp_arr_name()
         for d,nm in zip([mean_chem ,cell,mean_chg,rdf,adfa,adfb,ddf,nn],['mean_chem' ,'cell','mean_chg','rdf','adfa','adfb','ddf','nn']):
          if d!=[]:
            for ff,dd in enumerate(d):
             cat.append(dd)
             if nm=='mean_chem':
                tag=chem_nm[ff]
             else:
               tag=str(nm)+str('_')+str(ff)
             nmes.append(str(tag))
         cat=np.array(cat).astype(float) 
         #print (nmes,len(nmes))
         return nmes
        else:
         for d,nm in zip([mean_chem ,cell,mean_chg,rdf,adfa,adfb,ddf,nn],['mean_chem' ,'cell','mean_chg','rdf','adfa','adfb','ddf','nn']):
          if d!=[]:
            for ff,dd in enumerate(d):
             cat.append(dd)
         cat=np.array(cat).astype(float) 
      ##except:
       # pass 
        return cat

def get_chemonly(string=''):
  comp=Composition(string)
  el_dict=comp.get_el_amt_dict()
  arr=[]
  for k,v in el_dict.items():
     des=get_descrp_arr(k)
     arr.append(des)
  mean_chem=np.mean(arr,axis=0)
  return mean_chem

if __name__ == '__main__':
  
  #chemo-structural features
  s=Structure.from_file('POSCAR')
  x=get_comp_descp(s)
  print (len(x))

  #only chemical features for a structure
  y=get_comp_descp(struct=s,jcell=False,jmean_chem=True,jmean_chg=False,jrdf=False,jrdf_adf=False,print_names=False)
  print (len(y))


  #chemical features based on composition only
  chm=get_chemonly('Al2O3')
  print (len(chm))
   
  #charge descriptors
  comp=Composition('Al2O3')
  el_dict=comp.get_el_amt_dict()
  chgarr=[]
  for k,v in el_dict.items():
     chg=get_chgdescrp_arr(k)
     chgarr.append(chg)
  mean_chg=np.mean(chgarr,axis=0)
  print (len(mean_chg))

