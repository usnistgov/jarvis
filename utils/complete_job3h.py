import sys

import numpy as np
from ase.calculators.vasp import VaspChargeDensity
from monty.serialization import loadfn, dumpfn
from monty.json import MontyDecoder,MontyEncoder
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Vasprun
import sys,json,os,glob
from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
#provide your own MAPI Key from pymatgen
#MPall_datacopy.json
MAPI_KEY = "xZwPL7aKBHq3Gehb" 
info=[]
def getprops(dir=''):
    pwd=os.getcwd()
    os.chdir(dir)
    mbandgap='na'
    bbandgap='na'
    c11='na'
    n300='na'
    phi='na'
    p300='na'
    unit='Gpa'
    scl=1.0
    for sb_file in glob.glob("*.json"):
        if "MAIN-RELAX" in sb_file:
            #jfile=str(os.getcwd())+str("/")+str(sb_file)
            try:
              m_xml=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/vasprun.xml")
              mkpfile=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/KPOINTS")
   
              v_xml=Vasprun(m_xml)
              bands = v_xml.get_band_structure(mkpfile, line_mode =False, efermi = v_xml.efermi)
              if bands.get_band_gap()['direct']==True:
                    mbandgap=str(round(bands.get_band_gap()['energy'],3))+str(" D")
              else:
                    mbandgap=str(round(bands.get_band_gap()['energy'],3))+str(" I")


              LOCPOTfile=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/LOCPOT")
              allowed = "xyzXYZ"
              direction = 'z'#sys.argv[2].lstrip()
              #if allowed.find(direction) == -1 or len(direction)!=1 :
              #    print "** WARNING: The direction was input incorrectly."
              #    print "** Setting to z-direction by default."
              if direction.islower():
                  direction = direction.upper()
              filesuffix = "_%s" % direction

              # Open geometry and density class objects
              #-----------------------------------------
              vasp_charge = VaspChargeDensity(filename = LOCPOTfile)
              potl = vasp_charge.chg[-1]
              atoms = vasp_charge.atoms[-1]
              del vasp_charge

              # For LOCPOT files we multiply by the volume to get back to eV
              if 'LOCPOT' in LOCPOTfile:
                  potl=potl*atoms.get_volume()

              #print "\nReading file: %s" % LOCPOTfile
              #print "Performing average in %s direction" % direction

              # Read in lattice parameters and scale factor
              #---------------------------------------------
              cell = atoms.cell

              # Find length of lattice vectors
              #--------------------------------
              latticelength = np.dot(cell, cell.T).diagonal()
              latticelength = latticelength**0.5

              # Read in potential data
              #------------------------
              ngridpts = np.array(potl.shape)
              totgridpts = ngridpts.prod()
              #print "Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2])
              #print "Total number of points is %d" % totgridpts
              #print "Reading potential data from file...",
              sys.stdout.flush()
              #print "done."

              # Perform average
              #-----------------
              if direction=="X":
                  idir = 0
                  a = 1
                  b = 2
              elif direction=="Y":
                  a = 0
                  idir = 1
                  b = 2
              else:
                  a = 0
                  b = 1
                  idir = 2
              a = (idir+1)%3
              b = (idir+2)%3
              # At each point, sum over other two indices
              average = np.zeros(ngridpts[idir],np.float)
              for ipt in range(ngridpts[idir]):
                  if direction=="X":
                      average[ipt] = potl[ipt,:,:].sum()
                  elif direction=="Y":
                      average[ipt] = potl[:,ipt,:].sum()
                  else:
                      average[ipt] = potl[:,:,ipt].sum()

              if 'LOCPOT' in LOCPOTfile:
                  # Scale by number of grid points in the plane.
                  # The resulting unit will be eV.
                  average /= ngridpts[a]*ngridpts[b]
              else:
                  # Scale by size of area element in the plane,
                  # gives unit e/Ang. I.e. integrating the resulting
                  # CHG_dir file should give the total charge.
                  area = np.linalg.det([ (cell[a,a], cell[a,b] ),
                                         (cell[b,a], cell[b,b])])
                  dA = area/(ngridpts[a]*ngridpts[b])
                  average *= dA

              # Print out average
              #-------------------
              averagefile = LOCPOTfile + filesuffix
              #print "Writing averaged data to file %s..." % averagefile,
              sys.stdout.flush()
              outputfile = open(averagefile,"w")
              if 'LOCPOT' in LOCPOTfile:
                  outputfile.write("#  Distance(Ang)     Potential(eV)\n")
              else:
                  outputfile.write("#  Distance(Ang)     Chg. density (e/Ang)\n")
              xdiff = latticelength[idir]/float(ngridpts[idir]-1)
              xs=[]
              ys=[]
              for i in range(ngridpts[idir]):
                  x = i*xdiff
                  xs.append(x)
                  ys.append(average[i])
                  outputfile.write("%15.8g %15.8g\n" % (x,average[i]))
              outputfile.close()
              #print "done."
              #run=str(os.getcwd())+str("/")+str(fold.split(".json")[0])+str("/")+str("vasprun.xml")
              # run=str(os.getcwd())+str("/")+str(fold)+str("/vasprun.xml")
              vrun=v_xml#Vasprun(run)
              Ef=vrun.efermi
              avg_max=max(average)
              phi=float(avg_max)-float(Ef)


              #complete_dos = v_xml.complete_dos
              #mbandgap=float(complete_dos.get_gap())
            except:
               pass
        if "MAIN-BAND" in sb_file:
            #jfile=str(os.getcwd())+str("/")+str(sb_file)
            try:
              m_xml=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/vasprun.xml")
              bkpfile=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/KPOINTS")
              v_xml=Vasprun(m_xml)
              bands = v_xml.get_band_structure(bkpfile, line_mode =False, efermi = v_xml.efermi)
              if bands.get_band_gap()['direct']==True:
                    bbandgap=str(round(bands.get_band_gap()['energy'],3))+str(" D")
              else:
                    bbandgap=str(round(bands.get_band_gap()['energy'],3))+str(" I")
   
              #complete_dos = v_xml.complete_dos
              #bbandgap=float(complete_dos.get_gap())
            except:
               pass
        if "MAIN-ELASTIC" in sb_file:
            try:
               m_outcar=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/OUTCAR")
               outcar=Outcar(m_outcar)
               elast=(outcar.elastic_tensor)/float(10)
               if '1L' in dir:
                   sfile=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/POSCAR")
                   struct=Structure.from_file(sfile)
                   scl=0.1*float(abs(struct.lattice.matrix[2][2]))
                   unit='Nm'
               c11=scl*elast[0][0]
               
            except:
                pass
        if "Effective_Mass" in sb_file:
            try:
              efile=str(os.getcwd())+str("/")+str('Effective_Mass.json')
              eff= loadfn(efile, cls=MontyDecoder)
              n300=eff['n']['300'][0][0]
              p300=eff['p']['300'][0][0]
            except:
               pass
    os.chdir(pwd)
    c11=str(c11)+str(unit)
    #print "dir,phi=",dir,phi
    #print "bandgap=",mbandgap,bbandgap,c11,n300,p300
    return mbandgap,bbandgap,c11,n300,p300,phi

#bandgap,bbandgap,c11,n300,p300=getprops("/data/knc6/NIST4/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO")

#sys.exit()
def get_struct_from_mp():
    ratio=0
    mpids=[]
    pf='na'
    compat=False
    cnt=0
    #with open('MPall_datacopy.json') as json_data:
    #    data = json.load(json_data)i
    data = loadfn('/data/knc6/OQMD/MPall/MPall_datacopy.json', cls=MontyDecoder)
    for d in data:
            mpid= str(d['mp_id'])
            ini= (d['ini_structure'])
            fin= (d['fin_structure'])
            finder = SpacegroupAnalyzer(ini)
            latt=finder.get_lattice_type()
            bg=float(d["bg"])
            a,b,c=ini.lattice.abc
            ehull=d["ehull"]
            icsd=d["icsd"]
            sp=len(ini.types_of_specie)
            cnt=cnt+1
            a1,b1,c1=fin.lattice.abc
            ratioc= round(abs((float(c)-float(c1))/float(c1)),2)
            ratiob= round(abs((float(b)-float(b1))/float(b1)),2)
            ratioa= round(abs((float(a)-float(a1))/float(a1)),2)
            if  ehull==0 and latt!='cubic' and sp>1 and ratioc>0.01 and icsd!=None and icsd!=[] :
               info.append(mpid)  
    return info
def give_structure(id=''):
    data = loadfn('/data/knc6/OQMD/MPall/MPall_datacopy.json', cls=MontyDecoder)
    for d in data:
            mpid= str(d['mp_id'])
            ini= (d['ini_structure'])
            fin= (d['fin_structure'])
            if mpid==id:
                struct=fin
    return struct
def get_formula(mpid):
    ratio='na'
    rf='na'
    sg='na'
    af='na'
    ehull='na'
    a='na'
    b='na'
    c='na'
    a1='na'
    b1='na'
    c1='na'
    compat=False
    structure='na'
    ic_structure='na'
    try:
      with MPRester(MAPI_KEY) as m:
        data = m.get_data(mpid)
        cnt=0
      
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            bg=float(d["band_gap"])
            structure = m.get_structure_by_material_id(x['material_id'])
            af=str(structure.composition.anonymized_formula)
            rf=str(structure.composition.reduced_formula)
            a,b,c=structure.lattice.abc
            x['spacegroup'] = d['spacegroup']
            sg = x['spacegroup']['symbol']
            #sg = sg.replace('/', 'slash')
#           print  x['spacegroup'],x['material_id']
            cnt=cnt+1
            ic_structure = m.get_structure_by_material_id(x['material_id'],final=False)
            a1,b1,c1=ic_structure.lattice.abc
            ratio= round(((float(c)-float(c1))/float(c1)),2)
            ehull=d["e_above_hull"]
    except:
      pass
    return af,ehull,ratio,rf,sg,structure,ic_structure
    
import os,glob,sys
all_folds=[]
mpidlist=[]
enmpidlist=[]
bmpidlist=[]
enbmpidlist=[]
elmpidlist=[]
enelmpidlist=[]
smpidlist=[]
ensmpidlist=[]
TD1Lids=[]
TDids=[]
SSids=[]
EBids=[]
for fil in glob.glob("*-1L/*_*"):
    
#for fil in glob.glob("*-*/*_*"):
    id=str(fil).split('-')[2]+str('-')+str(fil).split('-')[3]
    if id not in TD1Lids:
       TD1Lids.append(id)
    if os.path.isdir(fil) :
        tag1=0
        tag2=0
        pwd=os.getcwd()
        dirp=str(pwd)+str("/")+str(fil)
        os.chdir(dirp)
        for file in glob.glob("*"):
            tag1=0
            tag2=0
            for sub_file in glob.glob("*"):
                 if "ENCUT" in sub_file  :
                    tag2=2
                 if "MAIN-BAND" in sub_file and 'PBEBO' in fil:
                   tag1=1
                   for sb_file in glob.glob("*.json"):
                       if "MAIN-RELAX" in sb_file:
                           case={}
                           jfile=str(os.getcwd())+str("/")+str(sb_file)
                           jdata = loadfn(jfile, cls=MontyDecoder)
                           mdir=str(os.getcwd())
                           case['dir']=mdir
                           case['mpid']=id
                           for data in jdata:
#                              #print data["final_energy"],data["poscar_final"]
                              case['enp']=float(data["final_energy"])/float(len(data["poscar_final"]))
                              case['poscar']=(data["poscar_final"])
                              if case not in enmpidlist:
                                  enmpidlist.append(case)

        if tag1==1 and tag2==2:
                  name=fil
                  #mpidlist.append(name)
                  mpidlist.append(id)
        os.chdir(pwd)


for fil in glob.glob("Elements-bulk/*_*"):
#for fil in glob.glob("*-*/*_*"):
    id=str(fil.split('/')[1]).split('_')[0]
    if id not in EBids:
       EBids.append(id)
    if os.path.isdir(fil) :
        tag1=0
        tag2=0
        pwd=os.getcwd()
        dirp=str(pwd)+str("/")+str(fil)
        os.chdir(dirp)
        for file in glob.glob("*"):
            tag1=0
            tag2=0
            for sub_file in glob.glob("*"):
                 if "ENCUT" in sub_file  :
                    tag2=2
                 if "MAIN-BAND" in sub_file and 'PBEBO' in fil:
                   tag1=1
                   for sb_file in glob.glob("*.json"):
                       if "MAIN-RELAX" in sb_file:
                           case={}
                           jfile=str(os.getcwd())+str("/")+str(sb_file)
                           jdata = loadfn(jfile, cls=MontyDecoder)
                           case['mpid']=id
                           mdir=str(os.getcwd())
                           case['dir']=mdir
                           for data in jdata:
#                              #print data["final_energy"],data["poscar_final"]
                              case['enp']=float(data["final_energy"])/float(len(data["poscar_final"]))
                              case['poscar']=(data["poscar_final"])
                              if case not in enelmpidlist:
                                  enelmpidlist.append(case)

        if tag1==1 and tag2==2:
                  name=fil
                  #mpidlist.append(name)
                  elmpidlist.append(id)
        os.chdir(pwd)




for fil in glob.glob("2D-bulk/*_*"):
#for fil in glob.glob("*-*/*_*"):
    id=str(fil.split('/')[1]).split('_')[0]
    if id not in TDids:
       TDids.append(id)
    if os.path.isdir(fil) :
        tag1=0
        tag2=0
        pwd=os.getcwd()
        dirp=str(pwd)+str("/")+str(fil)
        os.chdir(dirp)
        for file in glob.glob("*"):
            tag1=0
            tag2=0
            for sub_file in glob.glob("*"):
                 if "ENCUT" in sub_file  :
                    tag2=2
                 if "MAIN-BAND" in sub_file and 'PBEBO' in fil:
                   tag1=1
                   for sb_file in glob.glob("*.json"):
                       if "MAIN-RELAX" in sb_file:
                           case={}
                           jfile=str(os.getcwd())+str("/")+str(sb_file)
                           jdata = loadfn(jfile, cls=MontyDecoder)
                           case['mpid']=id
                           mdir=str(os.getcwd())
                           case['dir']=mdir
                           for data in jdata:
#                              print data["final_energy"],data["poscar_final"]
                              case['enp']=float(data["final_energy"])/float(len(data["poscar_final"]))
                              case['poscar']=(data["poscar_final"])
                              if case not in enbmpidlist:
                                  enbmpidlist.append(case)

        if tag1==1 and tag2==2:
                  name=fil
                  #mpidlist.append(name)
                  bmpidlist.append(id)
        os.chdir(pwd)
for fil in glob.glob("Solar-*/*_*"):
#for fil in glob.glob("*-*/*_*"):
    id=str(fil.split('/')[1]).split('_')[0]
    if id not in SSids:
        SSids.append(id)
    if os.path.isdir(fil) :
        tag1=0
        tag2=0
        pwd=os.getcwd()
        dirp=str(pwd)+str("/")+str(fil)
        os.chdir(dirp)
        for file in glob.glob("*"):
            tag1=0
            tag2=0
            for sub_file in glob.glob("*"):
                 if "ENCUT" in sub_file  :
                    tag2=2
                 if "MAIN-BAND" in sub_file and 'PBEBO' in fil:
                   tag1=1
                   for sb_file in glob.glob("*.json"):
                       if "MAIN-RELAX" in sb_file:
                           case={}
                           jfile=str(os.getcwd())+str("/")+str(sb_file)
                           jdata = loadfn(jfile, cls=MontyDecoder)
                           case['mpid']=id
                           mdir=str(os.getcwd())
                           case['dir']=mdir
                           for data in jdata:
#                              #print data["final_energy"],data["poscar_final"]
                              case['enp']=float(data["final_energy"])/float(len(data["poscar_final"]))
                              case['poscar']=(data["poscar_final"])
                              if case not in ensmpidlist:
                                 ensmpidlist.append(case)

        if tag1==1 and tag2==2:
                  name=fil
                  #mpidlist.append(name)
                  smpidlist.append(id)
        os.chdir(pwd)

print "From 2D-1L", mpidlist,len(mpidlist)
print "From 2D-bulk",bmpidlist,len(bmpidlist)
print "From Solar-Semi",smpidlist,len(smpidlist)
arr=[]
dumpa=[]
diffallcount=0
diffcount=0
criteria=[]
for i in mpidlist:
    if i in bmpidlist :
       case={}
       arr.append(i)
       for a in enmpidlist:
           if a['mpid']==i:
              enlayer=a['enp']
              poslayer=a['poscar']
              ldir=a['dir']
       for b in enbmpidlist:
           if b['mpid']==i:
              enbulk=b['enp']
              posbulk=b['poscar']
              bdir=b['dir']
       diff=enlayer-enbulk
       case['mpid']=i
       case['diff']=diff

       diffallcount=diffallcount+1
       if float(diff)<0.2:
          diffcount=diffcount+1

       af,ehull,ratio,rf,sg,mpstr,icstr=get_formula(i)
       case['af']=af
       case['rf']=rf
       case['sg']=str(sg)
       case['ehull']=ehull
       case['ratio']=ratio
       case['poslayer']=poslayer
       case['mpstr']=mpstr
       case['icstr']=icstr
       case['posbulk']=posbulk
       case['ldir']=ldir
       case['bdir']=bdir
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(ldir)
       case['lmbandgap']=mbandgap
       case['lphi']=phi
       case['lbbandgap']=bbandgap
       case['lc11']=c11
       case['ln300']=n300
       case['lp300']=p300
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(bdir)
       case['bphi']=phi
       case['bmbandgap']=mbandgap
       case['bbbandgap']=bbandgap
       case['bc11']=c11
       case['bn300']=n300
       case['bp300']=p300
       if case not in criteria:
          criteria.append(case)
    elif i in smpidlist :
       case={}
       arr.append(i)
       for a in enmpidlist:
           if a['mpid']==i:
              enlayer=a['enp']
              poslayer=a['poscar']
              ldir=a['dir']
       for b in ensmpidlist:
           if b['mpid']==i:
              enbulk=b['enp']
              posbulk=b['poscar']
              bdir=b['dir']
       diff=enlayer-enbulk
       diffallcount=diffallcount+1
       if float(diff)<0.2:
          diffcount=diffcount+1
       case['mpid']=i
       case['diff']=diff
       af,ehull,ratio,rf,sg,mpstr,icstr=get_formula(i)
       case['mpstr']=mpstr
       case['icstr']=icstr
       case['af']=af
       case['rf']=rf
       case['sg']=str(sg)
       case['ehull']=ehull
       case['ratio']=ratio
       case['poslayer']=poslayer
       case['posbulk']=posbulk
       case['ldir']=ldir
       case['bdir']=bdir
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(ldir)
       case['lphi']=phi
       case['lmbandgap']=mbandgap
       case['lbbandgap']=bbandgap
       case['lc11']=c11
       case['ln300']=n300
       case['lp300']=p300
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(bdir)
       case['bphi']=phi
       case['bmbandgap']=mbandgap
       case['bbbandgap']=bbandgap
       case['bc11']=c11
       case['bn300']=n300
       case['bp300']=p300
       if case not in criteria:
          criteria.append(case)
    elif i in elmpidlist :
       case={}
       arr.append(i)
       for a in enmpidlist:
           if a['mpid']==i:
              enlayer=a['enp']
              poslayer=a['poscar']
              ldir=a['dir']
       for b in enelmpidlist:
           if b['mpid']==i:
              enbulk=b['enp']
              posbulk=b['poscar']
              bdir=b['dir']
       diff=enlayer-enbulk
       diffallcount=diffallcount+1
       if float(diff)<0.2:
          diffcount=diffcount+1
       case['mpid']=i
       case['diff']=diff
       af,ehull,ratio,rf,sg,mpstr,icstr=get_formula(i)
       case['mpstr']=mpstr
       case['icstr']=icstr
       case['af']=af
       case['rf']=rf
       case['sg']=str(sg)
       case['ehull']=ehull
       case['ratio']=ratio
       case['poslayer']=poslayer
       case['posbulk']=posbulk
       case['ldir']=ldir
       case['bdir']=bdir
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(ldir)
       case['lphi']=phi
       case['lmbandgap']=mbandgap
       case['lbbandgap']=bbandgap
       case['lc11']=c11
       case['ln300']=n300
       case['lp300']=p300
       mbandgap,bbandgap,c11,n300,p300,phi=getprops(bdir)
       case['bphi']=phi
       case['bmbandgap']=mbandgap
       case['bbbandgap']=bbandgap
       case['bc11']=c11
       case['bn300']=n300
       case['bp300']=p300
       if case not in criteria:
          criteria.append(case)
    else:
       dumpa.append(i)

print "COMMON BULK and Single Layer",arr,len(arr)
print
print
"""
info=get_struct_from_mp()
remain1=[]
for i in  bmpidlist:
    if i in info and i  not in arr:
       remain1.append(i)
print "remain1 try for 2D 1L",remain1,len(remain1)
print
print
remain2=[]
for i in  smpidlist:
    if i in info and i  not in arr:
       remain2.append(i)
"""
#f=open('Criteria3f.json','w')
f=open('Criteria3g.json','w')
f.write(json.dumps(criteria,indent=4,cls=MontyEncoder ))
f.close()
#print "remain2 from Solar-semi try for 2D 1L",remain2,len(remain2)
print
print
#print "criteria",criteria
print
print
print "enmp",enmpidlist
print
print
print "enbmp",enbmpidlist
print
print
print "ensmp",ensmpidlist
print
print "TD1Lids",TD1Lids,len(TD1Lids)
print
print "TDids",TDids,len(TDids)
print
print "SSids",SSids,len(SSids)
print
print "diffallcount",diffallcount
print
print "diffcount",diffcount
