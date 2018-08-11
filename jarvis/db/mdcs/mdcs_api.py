# Make an account at https://jarvis.nist.gov/
# MDCS Github: https://github.com/faical-yannick-congo/MDCS, https://github.com/ztrautt/MDCS
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.serialization import loadfn, dumpfn
import os,glob,json,mdcs,sys
from monty.json import MontyEncoder, MontyDecoder
import numpy as np
from pymatgen.io.vasp.outputs import Vasprun
from mdcs import curate,explore
from pymatgen.core.structure import Structure
#import dicttoxml
from xml.etree.ElementTree import Element, SubElement, Comment, tostring, ElementTree

f=open('passwrd.txt','r') #path to your passowrd file
passd=f.read().splitlines()[0]
f.close()

user='abc' #your username
schema='5a7f2a872887d500ab0c0d02'

# Converting data to xml
def data_json(energy='na',typ='na',formula='na',sgp='na',name='na',ref='na',func='na',elem='na',encut='na',kpoints='na',el_tens='na',KV='na',GV='na',m_eg='na',b_eg='na',op_eg='na',mbj_eg='na',en_arr='na',realx_arr='na',imagx_arr='na',realy_arr='na',imagy_arr='na',realz_arr='na',imagz_arr='na',men_arr='na',mrealx_arr='na',mimagx_arr='na',mrealy_arr='na',mimagy_arr='na',mrealz_arr='na',mimagz_arr='na',struct='na',other='na'):
    
    top = Element('JARVIS-DFT')
    child = SubElement(top, 'energy-ev')
    child.text = str(energy)
    child = SubElement(top, 'formula')
    child.text = str(formula)
    child = SubElement(top, 'space-group')
    child.text = str(sgp)
    child = SubElement(top, 'JVID')
    child.text = str(name)
    child = SubElement(top, 'reference')
    child.text = str(ref)
    child = SubElement(top, 'calc_type')
    child.text = str(typ)
    child = SubElement(top, 'functional')
    child.text = str(func)
    child = SubElement(top, 'elements')
    child.text = str(elem)
    child = SubElement(top, 'encut-ev')
    child.text = str(encut)
    child = SubElement(top, 'kpoints')
    child.text = str(kpoints)
    child = SubElement(top, 'elastic_tensor-gpa')
    child.text = str(el_tens)
    child = SubElement(top, 'kv-gpa')
    child.text = str(KV)
    child = SubElement(top, 'gv-gpa')
    child.text = str(GV)
    child = SubElement(top, 'scf_eg-ev')
    child.text = str(m_eg)
    child = SubElement(top, 'brill_eg-ev')
    child.text = str(b_eg)
    child = SubElement(top, 'optics_eg-ev')
    child.text = str(op_eg)
    child = SubElement(top, 'mbj-eg_ev')
    child.text = str(mbj_eg)
    child = SubElement(top, 'opt_energy_arr-ev')
    child.text = str(en_arr)
    child = SubElement(top, 'realx_arr')
    child.text = str(realx_arr)
    child = SubElement(top, 'imagx_arr')
    child.text = str(imagx_arr)
    child = SubElement(top, 'realy_arr')
    child.text = str(realy_arr)
    child = SubElement(top, 'imagy_arr')
    child.text = str(imagy_arr)
    child = SubElement(top, 'realz_arr')
    child.text = str(realz_arr)
    child = SubElement(top, 'imagz_arr')
    child.text = str(imagz_arr)
    child = SubElement(top, 'mbj_energy_arr-ev')
    child.text = str(men_arr)
    child = SubElement(top, 'mbj_realx_arr')
    child.text = str(mrealx_arr)
    child = SubElement(top, 'mbj_imagx_arr')
    child.text = str(mimagx_arr)
    child = SubElement(top, 'mbj_realy_arr')
    child.text = str(mrealy_arr)
    child = SubElement(top, 'mbj_imagy_arr')
    child.text = str(mimagy_arr)
    child = SubElement(top, 'mbj_realz_arr')
    child.text = str(mrealz_arr)
    child = SubElement(top, 'mbj_imagz_arr')
    child.text = str(mimagz_arr)
    child = SubElement(top, 'structure')
    child.text = str(struct)
    child = SubElement(top, 'other')
    child.text = str(other)
    filename=str(name)+str('.xml')
    ElementTree(top).write(filename)
    curate(filename,filename,schema,'https://jarvis.nist.gov/',user,passd,cert=False)

# Get a particular record
def get_record(file=''):
   r=explore.select('https://jarvis.nist.gov/',user,passd,cert=False,title=file,format='json')
   return r


# Delete all record
def delete_all(file=''):
    r=explore.select_all('https://jarvis.nist.gov/',user,passd,cert=False,format='json')
    for i in r:
      id=i['_id']
      explore.delete(id,'https://jarvis.nist.gov/',user,passd,cert=False)



# Download 3D material data from https://figshare.com/articles/jdft_3d-7-7-2018_json/6815699
# Download 2D material data from https://figshare.com/articles/jdft_2d-7-7-2018_json/6815705 
# Curating JARVIS-DFT data 

d=loadfn('jdft_3d-7-7-2018.json',cls=MontyDecoder)


count=0
for i in d:
 filname=str(i['jid'])+str('.xml')
 if  not os.path.exists(filname) :
   count=count+1
   energy=str(i['fin_en'])+str(',')+str(i['form_enp'])
   formula=str(i['final_str'].composition.reduced_formula)
   sgp= str(SpacegroupAnalyzer(i['final_str']).get_spacegroup_symbol())
   name=str(i['jid'])
   print (name)
   ref=str(i['mpid'])
   func=str('OptB88vdW')
   elem=''
   species=list(set(i['final_str'].species))
   for j in species:
     elem=str(elem)+str(j.symbol)+str('-')
   encut=str(i['encut'])
   kpoints=str(i['kpoints'].kpts[0][0])+str('x')+str(i['kpoints'].kpts[0][1])+str('x')+str(i['kpoints'].kpts[0][2])
   el_tens=str(i['elastic'])
   KV=str(i['kv'])
   GV=str(i['gv'])
   op_eg=str(i['op_gap'])
   mbj_eg=str(i['mbj_gap'])
   realx_arr=str(i['epsx'])
   mrealx_arr=str(i['mepsx'])
   realy_arr=str(i['epsy'])
   mrealy_arr=str(i['mepsy'])
   realz_arr=str(i['epsz'])
   mrealz_arr=str(i['mepsz'])
   typ=str('3D')
   other=str('Citation: 1) DOI:10.1038/s41598-017-05402-0, 2) DOI: 10.1038/sdata.2018.82, 3) arXiv:1804.01033v2 ')
   struct=str(i['final_str'].to(fmt='cif').replace('# generated using pymatgen','# JARVIS-DFT')) 
   data_json(other=other,energy=energy,typ=typ,formula=formula,sgp=sgp,name=name,func=func,elem=elem,encut=encut,kpoints=kpoints,el_tens=el_tens,KV=KV,GV=GV,op_eg=op_eg,mbj_eg=mbj_eg,realx_arr=realx_arr,mrealx_arr=mrealx_arr,realy_arr=realy_arr,mrealy_arr=mrealy_arr,realz_arr=realz_arr,mrealz_arr=mrealz_arr,struct=struct)



#data_json()
#x=get_record(file='JVASP-48137.xml')[0]['content']['JARVIS-DFT']['structure']
print (Structure.from_str(x,fmt='cif'))
#delete_all()
#curate_data()
#delete_all()
#for k,v in r[0].items():
#for k,v in r[0]['content'].items():
#    print k
