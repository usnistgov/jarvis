from math import pi
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
def get_adf(s=None,cutoff=20.0,intvl=0.1):
    neighbors_lst = s.get_all_neighbors(cutoff)
    all_distances = np.concatenate(tuple(map(lambda x: \
     [itemgetter(1)(e) for e in x], neighbors_lst)))
    rdf_dict = {}
    dist_hist, dist_bins = np.histogram(all_distances, \
     bins=np.arange(0, cutoff + intvl, intvl), density=False)
    shell_vol = 4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) \
    - np.power(dist_bins[:-1], 3))
    number_density = s.num_sites / s.volume
    rdf = dist_hist / shell_vol / number_density
    return [{'distances': dist_bins[:-1], 'distribution': rdf}]

def VoroCoord(s=None):
    coord_finder = VoronoiCoordFinder(s)
    arr=[]
    for i,j in enumerate(s.sites):
       coord_no = int(round(coord_finder.get_coordination_number(i)))
       arr.append(coord_no)
    mean_cord=np.mean(np.array(arr))
    #print s.composition.formula,mean_cord
    #d=get_erdf(struct=s,val='max')
    #print 'erdf',d
    #d=get_erdf(struct=s,val='min')
    #print 'erdf',d
    #d=get_rdf(s=s)
    #print 'rdf',d
    return  mean_cord
def get_vor_coord(el='Na'):
    f=open('/rk2/knc6/JARVIS-DFT/MDCS/make_desciptors/uni.json','r')
    u=json.load(f,cls=MontyDecoder )
    f.close()
    c='na'
    key='MAIN-RELAX'
    for i,j in u.iteritems():
      #print i,j
       
      if str(i)==str(el):
          path=(j['file']) 
          for j in os.listdir( path ):
           if key in j and '.json' in j :
             #print path
             #tmp=path.replace('MAIN-RELAX',key)
             fil=str(path)+str("/")+str(j)
             s=loadfn(fil)[0]['contcar']
             #print s
             #s=Structure.from_file(str(path)+str("/")+str(j.split('.json')[0])+str("/CONTCAR"))
             c=VoroCoord(s)
             break
    return c

def get_gap(vasp_xml):

    tree = ET.ElementTree(file=vasp_xml)

    eigs = list(tree.iter(tag='eigenvalues'))[-1]
    for event, elem in ET.iterparse(vasp_xml):
      tag = elem.tag
      if tag=='dos':
        ef=float(elem.find("i").text)
    eig_sets = [eig_set.getchildren()
                for eig_set in eigs.iter('set')
                if ('comment' in eig_set.attrib
                    and 'kpoint' in eig_set.attrib['comment'])]

    def eig_set_to_lists(eig_set):
        return list(map(lambda el: list(map(float, el.text.split())), eig_set))

    eig_sets = list(map(eig_set_to_lists, eig_sets))
    all_eigs = [eig for kpt_eigs in eig_sets for eig in kpt_eigs]
    cb=[]
    vb=[]
    for i in all_eigs:
        l=i[0]
        if l>ef:
             cb.append(l)
        else:
             vb.append(l)
    gap=round(float(sorted(cb)[0]-sorted(vb)[-1]),3)
    return gap

def get_elast_tensor(path='',key='MAIN-ELAST'):
      el_tens='na'
      KV='na'
      GV='na'
      try:
        for j in os.listdir( path ):
          if key in j and '.json' in j :
             tmp=path.replace('MAIN-RELAX',key)
             out=str(tmp)+str("/")+str(j.split('.json')[0])+str("/OUTCAR")
        ratio_c=1.0
        v=open(out,"r")
        lines = v.read().splitlines()
        for i,line in enumerate(lines):
            if "TOTAL ELASTIC MODULI (kBar)" in  line:
                #nbands=int(line.split(">")[1].split("<")[0])
                c11= lines[i+3].split()[1]
                c12= lines[i+3].split()[2]
                c13= lines[i+3].split()[3]
                c14= lines[i+3].split()[4]
                c15= lines[i+3].split()[5]
                c16= lines[i+3].split()[6]
                c21= lines[i+4].split()[1]
                c22= lines[i+4].split()[2]
                c23= lines[i+4].split()[3]
                c24= lines[i+4].split()[4]
                c25= lines[i+4].split()[5]
                c26= lines[i+4].split()[6]
                c31= lines[i+5].split()[1]
                c32= lines[i+5].split()[2]
                c33= lines[i+5].split()[3]
                c34= lines[i+5].split()[4]
                c35= lines[i+5].split()[5]
                c36= lines[i+5].split()[6]
                c41= lines[i+6].split()[1]
                c42= lines[i+6].split()[2]
                c43= lines[i+6].split()[3]
                c44= lines[i+6].split()[4]
                c45= lines[i+6].split()[5]
                c46= lines[i+6].split()[6]
                c51= lines[i+7].split()[1]
                c52= lines[i+7].split()[2]
                c53= lines[i+7].split()[3]
                c54= lines[i+7].split()[4]
                c55= lines[i+7].split()[5]
                c56= lines[i+7].split()[6]
                c61= lines[i+8].split()[1]
                c62= lines[i+8].split()[2]
                c63= lines[i+8].split()[3]
                c64= lines[i+8].split()[4]
                c65= lines[i+8].split()[5]
                c66= lines[i+8].split()[6]
                c11=round(ratio_c*float(c11)/float(10),1)
                c12=round(ratio_c*float(c12)/float(10),1)
                c13=round(ratio_c*float(c13)/float(10),1)
                c14=round(ratio_c*float(c14)/float(10),1)
                c15=round(ratio_c*float(c15)/float(10),1)
                c16=round(ratio_c*float(c16)/float(10),1)
                c21=round(ratio_c*float(c21)/float(10),1)
                c22=round(ratio_c*float(c22)/float(10),1)
                c23=round(ratio_c*float(c23)/float(10),1)
                c24=round(ratio_c*float(c24)/float(10),1)
                c25=round(ratio_c*float(c25)/float(10),1)
                c26=round(ratio_c*float(c26)/float(10),1)
                c31=round(float(c31)/float(10),1)
                c32=round(float(c32)/float(10),1)
                c33=round(float(c33)/float(10),1)
                c34=round(float(c34)/float(10),1)
                c35=round(float(c35)/float(10),1)
                c36=round(float(c36)/float(10),1)
                c41=round(float(c41)/float(10),1)
                c42=round(float(c42)/float(10),1)
                c43=round(float(c43)/float(10),1)
                c44=round(float(c44)/float(10),1)
                c45=round(float(c45)/float(10),1)
                c46=round(float(c46)/float(10),1)
                c51=round(float(c51)/float(10),1)
                c52=round(float(c52)/float(10),1)
                c53=round(float(c53)/float(10),1)
                c54=round(float(c54)/float(10),1)
                c55=round(float(c55)/float(10),1)
                c56=round(float(c56)/float(10),1)
                c61=round(float(c61)/float(10),1)
                c62=round(float(c62)/float(10),1)
                c63=round(float(c63)/float(10),1)
                c64=round(float(c64)/float(10),1)
                c65=round(float(c65)/float(10),1)
                c66=round(float(c66)/float(10),1)
                KV=float((c11+c22+c33)+2*(c12+c23+c31))/float(9)
                GV=float((c11+c22+c33)-(c12+c23+c31)+3*(c44+c55+c66))/float(15)
                KV=str(round(KV,3))
                GV=str(round(GV,3))
                break
        v.close()
        el_tens=str(c11)+str(',')+str(c12)+str(',')+str(c13)+str(',')+str(c14)+str(',')+str(c15)+str(',')+str(c16)+str(',')+str(c21)+str(',')+str(c22)+str(',')+str(c23)+str(',')+str(c24)+str(',')+str(c25)+str(',')+str(c26)+str(',')+str(c31)+str(',')+str(c32)+str(',')+str(c33)+str(',')+str(c34)+str(',')+str(c35)+str(',')+str(c36)+str(',')+str(c41)+str(',')+str(c42)+str(',')+str(c43)+str(',')+str(c44)+str(',')+str(c45)+str(',')+str(c46)+str(',')+str(c51)+str(',')+str(c52)+str(',')+str(c53)+str(',')+str(c54)+str(',')+str(c55)+str(',')+str(c56)+str(',')+str(c61)+str(',')+str(c62)+str(',')+str(c63)+str(',')+str(c64)+str(',')+str(c65)+str(',')+str(c66)
      except:
         pass

      #print 'el_tens',np.array(el_tens.split(',')).astype(float)
      return el_tens,KV,GV


def optoelelectronic(path='',key='MAIN-MBJ'):
    cwd=str(os.getcwd())
    op_gap='na'
    realx_arr='na'
    realy_arr='na'
    realz_arr='na'
    en_arr='na'
    #print 'path1',path
    try:
      for j in os.listdir( path ):
          if key in j and '.json' in j :
             #print path
             #tmp=path.replace('MAIN-RELAX',key)
             vrun=str(path)+str("/")+str(j.split('.json')[0])+str("/vasprun.xml")
             #op_gap=get_gap(vrun)
             #print vrun,op_gap
      run = Vasprun(vrun)
         
      op_gap=run.eigenvalue_band_properties[0]
      #print 'path',path,op_gap
      erange=len(run.dielectric[0])

      en=[]
      realx=[]
      #imagx=[]
      realy=[]
      #imagy=[]
      realz=[]
      #imagz=[]
      for i in range(0,erange-1):
        #en.append(str(run.dielectric[0][i]))
        realx.append(str(run.dielectric[1][i][0]))
        #imagx.append(str(run.dielectric[2][i][0]))
        realy.append(str(run.dielectric[1][i][1]))
        #imagy.append(str(run.dielectric[2][i][1]))
        realz.append(str(run.dielectric[1][i][2]))
        #imagz.append(str(run.dielectric[2][i][2]))

      realx_arr=realx[0]#','.join(realx)
      #imagx_arr=imagx[0]#','.join(imagx)
      realy_arr=realy[0]#','.join(realy)
      #imagy_arr=imagy[0]#','.join(imagy)
      realz_arr=realz[0]#','.join(realz)
      #imagz_arr=imagz[0]#','.join(imagz)

      #print 'path',path,op_gap,realx_arr,realy_arr,realz_arr

    except:
       pass
    #print path,op_gap


    os.chdir(cwd)
    return op_gap,realx_arr,realy_arr,realz_arr



def coh_en(el=''):
  coh=open('cohesive_energies.json','r')
  coh_dat=json.load(coh)
  coh.close()
  return coh_dat[el]
def elec_eff(el=''):
  coh=open('elc_aff.json','r')
  coh_dat=json.load(coh,cls=MontyDecoder)
  coh.close()
  return coh_dat[el]

def uni_en_prop(el='Na'):
   en='na'
   mop_gap='na'
   mrealx_arr='na'
   mrealy_arr='na'
   mrealz_arr='na'
   op_gap='na'
   realx_arr='na'
   realy_arr='na'
   realz_arr='na'
   el_tens='na'
   coord='na'
   KV='na'
   GV='na'
   f=open('/rk2/knc6/JARVIS-DFT/MDCS/make_desciptors/uni.json','r')
   u=json.load(f,cls=MontyDecoder )
   f.close()
   for i,j in u.iteritems():
      #print i,j
      if str(i)==str(el):
          en=j['en']
          el_tens,KV,GV=get_elast_tensor(j['file']) # optoelelectronic(path='',key='MAIN-MBJ')
          op_gap,realx_arr,realy_arr,realz_arr=optoelelectronic(j['file'],key='MAIN-OPTICS')
          mop_gap,mrealx_arr,mrealy_arr,mrealz_arr=optoelelectronic(j['file'],key='MAIN-OPTICS')
          coord=get_vor_coord(el)
          #print str(el),j['file'],j['mpid'],el_tens,KV,GV
          #print op_gap,realx_arr,realy_arr,realz_arr
    #      if el=='N':
    #         print el,'found',en

   return en,op_gap,realx_arr,realy_arr,realz_arr,mop_gap,mrealx_arr,mrealy_arr,mrealz_arr,el_tens,KV,GV,coord
def uni_en(el='Na'):
   en='na'
   f=open('/rk2/knc6/JARVIS-DFT/MDCS/make_desciptors/uni.json','r')
   u=json.load(f,cls=MontyDecoder )
   f.close()
   for i,j in u.iteritems():
      #print i,j
      if str(i)==str(el):
          en=j['en']
          el_tens,KV,GV=get_elast_tensor(j['file'])
          op_gap,realx_arr,realy_arr,realz_arr=optoelelectronic(j['file'])
          #print str(el),j['file'],j['mpid'],el_tens,KV,GV
          print op_gap,realx_arr,realy_arr,realz_arr
    #      if el=='N':
    #         print el,'found',en

   return en

def log_to_numb(tag=''):
    if tag==True:
       return 1
    elif tag==False:
        return 0
    else:
        return 'na'
def block_num(tag=''):
    if tag=='s':
       return 1
    elif tag=='p':
       return 2
    elif tag=='d':
       return 3
    elif tag=='f':
       return 4
    else:
       return 0
def magfile_val(file='',atm_num=0):
     f=open(file,'r')
     lines=f.read().splitlines()
     f.close()
     val='na'
     for i,j in enumerate(lines):
         #print i,j,int(atm_num)
         if i+1==int(atm_num):
            val=j
            break
            #print 'printinh',i+1,atm_num,j.split()[0]
     #print 'val',val
     return val#j.split()[0]
def make_el_dat():
	el=Element
	#magfile_val('first_ion',1)
	#sys.exit()
	count=0
	el_data={}
	for i in el:
	 info={}
	 mem={}

	 #jv_enp,op_gap,realx_arr,realy_arr,realz_arr,mop_gap,mrealx_arr,mrealy_arr,\
         #  mrealz_arr,el_tens,KV,GV=uni_en_prop(str(i.symbol)) 

	 try:
	  #if i.symbol=='H':
	   
	  #if not i.is_actinoid:
	   count=count+1
          #if count<5:
	   zz=int(i.number)
	   symb=str(i.symbol)
	   ndun=magfile_val('NDunfilled',zz)
	   ndv=magfile_val('NDvalence',zz)
	   nsun=magfile_val('NSunfilled',zz)
	   nsv=magfile_val('NSvalence',zz)
	   npun=magfile_val('NPunfilled',zz)
	   npv=magfile_val('NPvalence',zz)
	   nfun=magfile_val('NFunfilled',zz)
	   nfv=magfile_val('NFvalence',zz)
	   ion=magfile_val('first_ion',zz)
	   oq_bg=magfile_val('oqmd_bg',zz)
	   el_ef=magfile_val('elec_aff',zz)
	   vol_pa=magfile_val('vol_pa',zz)
	   hfus=magfile_val('HFusion',zz)
	   oq_enp=magfile_val('oqmd_enp',zz)
	   polz=magfile_val('polariz',zz)
	   #jv_enp=uni_en(str(i.symbol)) 
	   

	   #mem['symbol']=str(i.symbol)
	   mem['Z']=int(zz)
	   mem['X']=float(i.X)
	   mem['ndunfill']=int(ndun)
	   mem['ndvalence']=int(ndv)
	   mem['nsunfill']=int(nsun)
	   mem['nsvalence']=int(nsv)
	   mem['npunfill']=int(npun)
	   mem['npvalence']=int(npv)
	   mem['nfunfill']=int(nfun)
	   mem['nfvalence']=int(nfv)
	   mem['first_ion_en']=float(ion)
	   mem['elec_aff']=float(el_ef)*0.01 #ev
	   mem['oq_bg']=float(oq_bg)
	   mem['oq_enp']=float(oq_enp)
	   mem['row']=int(i.row)
	   mem['coulmn']=int(i.group)
          
	   mem['max_oxid_s']=float(i.max_oxidation_state)
	   mem['min_oxid_s']=float(i.min_oxidation_state)
	   mem['block']=int(block_num(i.block))
	   mem['is_alkali']=int(log_to_numb(i.is_alkali))
	   mem['is_alkaline']=log_to_numb(i.is_alkaline)
	   mem['is_metalloid']=log_to_numb(i.is_metalloid)
	   mem['is_noble_gas']=log_to_numb(i.is_noble_gas)
	   mem['is_transition_metal']=log_to_numb(i.is_transition_metal)
	   mem['is_metalloid']=log_to_numb(i.is_metalloid)
	   mem['is_halogen']=log_to_numb(i.is_halogen)
	   mem['is_lanthanoid']=log_to_numb(i.is_lanthanoid)
	   mem['is_actinoid']=log_to_numb(i.is_actinoid)
	   mem['atom_mass']=float(str(i.atomic_mass).split('amu')[0])
	   mem['atom_rad']=float(str(i.atomic_radius).split('ang')[0])
	   mem['therm_cond']=float(str(i.thermal_conductivity).split('W m<sup>-1</sup> K<sup>-1</sup>')[0])
	   mem['mol_vol']=float(str(i.molar_volume).split('cm<sup>3</sup>')[0])
	   mem['bp']=float(str(i.boiling_point).split('K')[0])
	   mem['mp']=float(str(i.melting_point).split('K')[0])
	   mem['avg_ion_rad']=float(str(i.average_ionic_radius).split('ang')[0])
	   mem['hfus']=float(hfus)*float(0.01)
	   mem['polzbl']=float(polz)
           #print 'before'
	   jv_enp,op_gap,realx_arr,realy_arr,realz_arr,mop_gap,mrealx_arr,mrealy_arr,mrealz_arr,el_tens,KV,GV,coord=uni_en_prop(str(i.symbol)) 
           print coord
           #return en,op_gap,realx_arr,realy_arr,realz_arr,el_tens,KV,GV
	   #jv_enp=uni_en(str(i.symbol)) 
           cij=np.array(el_tens.split(',')).astype(float)
	   for ei,ej in enumerate(cij):
              el_coeff= str('C-')+str(ei)
              mem[el_coeff]=ej
           #print mem
           #import sys
           #sys.exit()
           mem['KV']=KV
           mem['voro_coord']=coord
           mem['GV']=GV
	   mem['jv_enp']=float(jv_enp)
           mem['op_eg']=op_gap 
           mem['mop_eg']=mop_gap 
           mem['e1']=realx_arr
           mem['e2']=realy_arr
           mem['e3']=realz_arr
           mem['me1']=mrealx_arr
           mem['me2']=mrealy_arr
           mem['me3']=mrealz_arr
	   #print ndun,nsv
           tmp=[]
           print 'iters1'
           iters=['hfus','polzbl','first_ion_en','elec_aff','mol_vol','bp','mp','atom_mass','atom_rad','therm_cond','X','voro_coord']
           for i1 in iters:
              for j1 in iters:
               tmp1=str(i1)+str('_add_')+str(j1)
               tmp2=str(j1)+str('_add_')+str(i1)
               if i1!=j1 and tmp1 not in tmp and tmp2 not in tmp :
                  #print i1,j1
                  kkey=str(i1)+str('_add_')+str(j1)
                  tmp.append(kkey)
                  vval=float(mem[i1])+float(mem[j1])
                  #print 'iters1',kkey
                  mem[kkey]=vval
		  kkey=str(i1)+str('_mult_')+str(j1)
		  vval=mem[i1]*mem[j1]
		  mem[kkey]=vval
                          
               if i1!=j1  :
		  kkey=str(i1)+str('_subs_')+str(j1)
		  vval=mem[i1]-mem[j1]
		  mem[kkey]=vval

           iters=['hfus','polzbl','first_ion_en','mol_vol','bp','mp','atom_mass','atom_rad','therm_cond','voro_coord']
           print 'iters2'
           for i1 in iters:
              for j1 in iters:
               if i1!=j1:
                  kkey=str(i1)+str('_divi_')+str(j1)
                  #print 'div key',kkey
                  vval=float(mem[i1])/float(mem[j1])
                  mem[kkey]=vval
             


	   #print "mem",mem,len(mem)
           #print
           #print
	   el_data[str(i.symbol)]=mem
           for k,v in mem.iteritems():
            if 'na' in v:
                print mem
	   #coh=coh_en(str(i.symbol)) 
	   #print i.symbol,i.number,i.row,i.group,i.Z,i.X, i.max_oxidation_state,i.min_oxidation_state,block_num(i.block),i.block,log_to_numb(i.is_alkali),log_to_numb(i.is_alkaline),log_to_numb(i.is_metalloid),log_to_numb(i.is_noble_gas),log_to_numb(i.is_transition_metal),log_to_numb(i.is_metalloid),log_to_numb(i.is_halogen),log_to_numb(i.is_lanthanoid),log_to_numb(i.is_actinoid),str(i.atomic_mass).split('amu')[0],str(i.atomic_radius).split('ang')[0],str(i.molar_volume).split('cm<sup>3</sup>')[0],str(i.boiling_point).split('K')[0],str(i.melting_point).split('K')[0],str(i.average_ionic_radius).split('ang')[0]
	   #print i.symbol,i.number,i.row,i.group,i.Z,i.X, i.max_oxidation_state,i.min_oxidation_state,block_num(i.block),i.block,log_to_numb(i.is_alkali),log_to_numb(i.is_alkaline),log_to_numb(i.is_metalloid),log_to_numb(i.is_noble_gas),log_to_numb(i.is_transition_metal),log_to_numb(i.is_metalloid),log_to_numb(i.is_halogen),log_to_numb(i.is_lanthanoid),log_to_numb(i.is_actinoid),str(i.atomic_mass).split('amu')[0],str(i.atomic_radius).split('ang')[0],str(i.molar_volume).split('cm<sup>3</sup>')[0],str(i.boiling_point).split('K')[0],str(i.melting_point).split('K')[0],str(i.average_ionic_radius).split('ang')[0],coh_dat[i.symbol]
	   #print i.symbol,i.row,i.group,i.Z,i.X, i.max_oxidation_state,i.min_oxidation_state,i.block,log_to_numb(i.is_alkali),log_to_numb(i.is_alkaline),log_to_numb(i.is_metalloid),log_to_numb(i.is_noble_gas),log_to_numb(i.is_transition_metal),log_to_numb(i.is_metalloid),log_to_numb(i.is_halogen),log_to_numb(i.is_lanthanoid),log_to_numb(i.is_actinoid),str(i.atomic_mass).split('amu')[0],str(i.atomic_radius).split('ang')[0],str(i.molar_volume).split('cm<sup>3</sup>')[0],str(i.boiling_point).split('K')[0],str(i.melting_point).split('K')[0],str(i.density_of_solid).split('kg m<sup>-3</sup>')[0],coh_dat[i.symbol]
	   #print i.symbol,i.row,i.group,i.Z,i.X, i.max_oxidation_state,i.min_oxidation_state,i.oxidation_states,i.block,i.is_alkali,i.is_alkaline,i.is_metalloid,i.is_noble_gas,i.is_transition_metal,i.is_metalloid,i.is_halogen,i.is_lanthanoid,i.is_actinoid,i.atomic_mass,i.atomic_radius,i.molar_volume,i.boiling_point,i.melting_point,i.density_of_solid,i.average_ionic_radius,i.ionic_radii
	 except:
	   pass
	f=open('Elements.json','w')
	f.write(json.dumps(el_data,indent=4))
	f.close() 
	print count
#make_el_dat()
def get_descrp_arr(elm=''):
      arr=[]
      try:      
	f=open('Elements.json','r')
	dat=json.load(f)
	f.close()
	d=dat[elm]
	arr=[]
	for k,v in d.iteritems():
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
    rdf = dist_hist / shell_vol / number_density
    return [{'distances': dist_bins[:-1], 'distribution': rdf}]

def get_erdf(struct=None,cutoff=20.0,dr=0.1,val='max'):
    struct = ValenceIonicRadiusEvaluator(struct).structure
    nbins = int(cutoff / dr) + 1
    redf_dict = {"distances": np.array([(i + 0.5) * dr for i in range(nbins)])\
    ,"distribution": np.zeros(nbins, dtype=np.float)}
    if val=='max':
     for site in struct.sites:
          this_charge = float(site.specie.max_oxidation_state)
          #this_charge = float(site.specie.oxi_state)
          neighs_dists = struct.get_neighbors(site, cutoff)
          for neigh, dist in neighs_dists:
             neigh_charge = float(neigh.specie.max_oxidation_state)
             #neigh_charge = float(neigh.specie.oxi_state)
             bin_index = int(dist / dr)
             redf_dict["distribution"][bin_index] += (this_charge * \
             neigh_charge) / (struct.num_sites * dist)
    else:
     for site in struct.sites:
          this_charge = float(site.specie.min_oxidation_state)
          #this_charge = float(site.specie.oxi_state)
          neighs_dists = struct.get_neighbors(site, cutoff)
          for neigh, dist in neighs_dists:
             neigh_charge = float(neigh.specie.min_oxidation_state)
             #neigh_charge = float(neigh.specie.oxi_state)
             bin_index = int(dist / dr)
             redf_dict["distribution"][bin_index] += (this_charge * \
             neigh_charge) / (struct.num_sites * dist)
    return [redf_dict]

def get_comp_descp(struct=''):  
       cat=[]
       try: 
        s=(struct).get_primitive_structure()
        a=s.lattice.abc[0]
        b=s.lattice.abc[1]
        c=s.lattice.abc[2]
        ab=float(s.lattice.abc[0])/float(s.lattice.abc[1])
        bc=float(s.lattice.abc[1])/float(s.lattice.abc[2])
        ca=float(s.lattice.abc[2])/float(s.lattice.abc[0])
        alph=s.lattice.angles[0]
        bet=s.lattice.angles[1]
        gam=s.lattice.angles[2]
        spg=str(SpacegroupAnalyzer(s).get_spacegroup_number()).split()[0]
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
	for k,v in el_dict.iteritems():
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
#make_el_dat()
#s=Structure.from_file('POSCAR')
#p=get_comp_descp(s)
#print p,len(p)
#get_erdf(struct=s,cutoff=20.0,dr=0.1)
#get_vor_coord(el='N')
#for i in Element:
  
#   get_vor_coord(i.symbol)
