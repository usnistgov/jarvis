import sys,math

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

def optics(ru=''):
    run = Vasprun(ru)
    erange=len(run.dielectric[0])
    en=[]
    realx=[]
    imagx=[]
    absorpx=[]
    refrx=[]
    reflx=[]
    eelsx=[]
    extcx=[]
    opt_conx=[]
    realy=[]
    imagy=[]
    absorpy=[]
    refry=[]
    refly=[]
    eelsy=[]
    extcy=[]
    opt_cony=[]
    realz=[]
    imagz=[]
    absorpz=[]
    refrz=[]
    reflz=[]
    eelsz=[]
    extcz=[]
    opt_conz=[]
    H=4.13566733*(10**(-15))
    #c0=2.99792458
    c0=2.99792458*(math.pow(10,8))
    for i in range(0,erange-1):
        en.append(run.dielectric[0][i])
        realx.append(run.dielectric[1][i][0])
        imagx.append(run.dielectric[2][i][0])

        ab_valx=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][0]+math.sqrt((run.dielectric[2][i][0])*(run.dielectric[2][i][0])+(run.dielectric[1][i][0])*(run.dielectric[1][i][0]))))/float(float(c0)*100.0)))
        #ab_valx=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][0]+math.sqrt((run.dielectric[2][i][0])*(run.dielectric[2][i][0])+(run.dielectric[1][i][0])*(run.dielectric[1][i][0]))))/float(float(c0)/100.0)))
        absorpx.append(ab_valx)
        refr_valx=float(math.sqrt((run.dielectric[1][i][0])+math.sqrt((run.dielectric[2][i][0])*(run.dielectric[2][i][0])+(run.dielectric[1][i][0])*(run.dielectric[1][i][0]))))/float(1.4142)
        refrx.append(refr_valx)
        eels_valx=float(run.dielectric[2][i][0])/float((run.dielectric[2][i][0])*(run.dielectric[2][i][0])+(run.dielectric[1][i][0])*(run.dielectric[1][i][0]))
        eelsx.append(eels_valx)
        extc_valx=float(math.sqrt(-(run.dielectric[1][i][0])+math.sqrt((run.dielectric[2][i][0])*(run.dielectric[2][i][0])+(run.dielectric[1][i][0])*(run.dielectric[1][i][0]))))/float(1.4142)
        extcx.append(extc_valx)
        refl_valx=((refr_valx-1)*(refr_valx-1)+extc_valx*extc_valx)/((refr_valx+1)*(refr_valx+1)+extc_valx*extc_valx)
        reflx.append(refl_valx)
        opt_valx=float(float(run.dielectric[0][i])/float(H))/float(4*math.pi)*(run.dielectric[2][i][0])
        opt_conx.append(opt_valx) #Real part of optical conductivity #http://www.wien2k.at/reg_user/textbooks/WIEN2k_lecture-notes_2013/optic_handout.pdf
        realy.append(run.dielectric[1][i][1])
        imagy.append(run.dielectric[2][i][1])

        ab_valy=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][1]+math.sqrt((run.dielectric[2][i][1])*(run.dielectric[2][i][1])+(run.dielectric[1][i][1])*(run.dielectric[1][i][1]))))/float(float(c0)*100.0)))
        #ab_valy=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][1]+math.sqrt((run.dielectric[2][i][1])*(run.dielectric[2][i][1])+(run.dielectric[1][i][1])*(run.dielectric[1][i][1]))))/float(float(c0)/100.0)))
        absorpy.append(ab_valy)
        refr_valy=float(math.sqrt((run.dielectric[1][i][1])+math.sqrt((run.dielectric[2][i][1])*(run.dielectric[2][i][1])+(run.dielectric[1][i][1])*(run.dielectric[1][i][1]))))/float(1.4142)
        refry.append(refr_valy)
        eels_valy=float(run.dielectric[2][i][1])/float((run.dielectric[2][i][1])*(run.dielectric[2][i][1])+(run.dielectric[1][i][1])*(run.dielectric[1][i][1]))
        eelsy.append(eels_valy)
        extc_valy=float(math.sqrt(-(run.dielectric[1][i][1])+math.sqrt((run.dielectric[2][i][1])*(run.dielectric[2][i][1])+(run.dielectric[1][i][1])*(run.dielectric[1][i][1]))))/float(1.4142)
        extcy.append(extc_valy)
        refl_valy=((refr_valy-1)*(refr_valy-1)+extc_valy*extc_valy)/((refr_valy+1)*(refr_valy+1)+extc_valy*extc_valy)
        refly.append(refl_valy)
        opt_valy=float(float(run.dielectric[0][i])/float(H))/float(4*math.pi)*(run.dielectric[2][i][1])
        opt_cony.append(opt_valy) #Real part of optical conductivity #http://www.wien2k.at/reg_user/textbooks/WIEN2k_lecture-notes_2013/optic_handout.pdf





        realz.append(run.dielectric[1][i][2])
        imagz.append(run.dielectric[2][i][2])

        ab_valz=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][2]+math.sqrt((run.dielectric[2][i][2])*(run.dielectric[2][i][2])+(run.dielectric[1][i][2])*(run.dielectric[1][i][2]))))/float(float(c0)*100.0)))
        #ab_valz=1.4142*((float(run.dielectric[0][i])/float(H))*(float(math.sqrt(-run.dielectric[1][i][2]+math.sqrt((run.dielectric[2][i][2])*(run.dielectric[2][i][2])+(run.dielectric[1][i][2])*(run.dielectric[1][i][2]))))/float(float(c0)/100.0)))
        absorpz.append(ab_valz)
        refr_valz=float(math.sqrt((run.dielectric[1][i][2])+math.sqrt((run.dielectric[2][i][2])*(run.dielectric[2][i][2])+(run.dielectric[1][i][2])*(run.dielectric[1][i][2]))))/float(1.4142)
        refrz.append(refr_valz)
        eels_valz=float(run.dielectric[2][i][2])/float((run.dielectric[2][i][2])*(run.dielectric[2][i][2])+(run.dielectric[1][i][2])*(run.dielectric[1][i][2]))
        eelsz.append(eels_valz)
        extc_valz=float(math.sqrt(-(run.dielectric[1][i][2])+math.sqrt((run.dielectric[2][i][2])*(run.dielectric[2][i][2])+(run.dielectric[1][i][2])*(run.dielectric[1][i][2]))))/float(1.4142)
        extcz.append(extc_valz)
        refl_valz=((refr_valz-1)*(refr_valz-1)+extc_valz*extc_valz)/((refr_valz+1)*(refr_valz+1)+extc_valz*extc_valz)
        reflz.append(refl_valz)
        opt_valz=float(float(run.dielectric[0][i])/float(H))/float(4*math.pi)*(run.dielectric[2][i][2])
        opt_conz.append(opt_valz) #Real part of optical conductivity #http://www.wien2k.at/reg_user/textbooks/WIEN2k_lecture-notes_2013/optic_handout.pdf

    opt_info={}
    try:
    		opt_info['en']=en[0]
    		opt_info['realx']=realx[0]
    		opt_info['imagx']=imagx[0]
    		opt_info['absorpx']=absorpx[0]
    		opt_info['refrx']=refrx[0]
    		opt_info['reflx']=reflx[0]
    		opt_info['eelsx']=eelsx[0]
    		opt_info['extcx']=extcx[0]
    		opt_info['opt_conx']=opt_conx[0]
    		opt_info['realy']=realy[0]
    		opt_info['imagy']=imagy[0]
    		opt_info['absorpy']=absorpy[0]
    		opt_info['refry']=refry[0]
    		opt_info['refly']=refly[0]
    		opt_info['eelsy']=eelsy[0]
    		opt_info['extcy']=extcy[0]
    		opt_info['opt_cony']=opt_cony[0]
    		opt_info['realz']=realz[0]
    		opt_info['imagz']=imagz[0]
    		opt_info['absorpz']=absorpz[0]
    		opt_info['refrz']=refrz[0]
    		opt_info['reflz']=reflz[0]
    		opt_info['eelsz']=eelsz[0]
    		opt_info['extcz']=extcz[0]
                opt_info['opt_conz']=opt_conz[0]
    except:
         pass
    return opt_info

def locpot(LOCPOTfile='',xml=None,direction='x'):
#def locpot(LOCPOTfile)='',direction='z'):
    #direction = 'z'#sys.argv[2].lstrip()
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
    vrun=xml#Vasprun(run)
    Ef=vrun.efermi
    avg_max=max(average)
    phi=float(avg_max)-float(Ef)
    return phi
info=[]


def getprops(dir=''):
    try:
      pwd=os.getcwd()
      os.chdir(dir)
      all_info={}
      mbandgap='na'
      bbandgap='na'
      name='na'
      c11='na'
      elast='na'
      eff='na'
      n300='na'
      phi='na'
      p300='na'
      unit='Gpa'
      scl=1.0
      for sb_file in glob.glob("*.json"):
        if "MAIN-OPTICS" in sb_file:
            try:
              o_xml=str(os.getcwd())+str("/")+str(str(sb_file).split(".json")[0])+str("/vasprun.xml")
              opt_info={}
              opt_info=optics(o_xml)
            except:
                pass
            all_info['opt_info']=opt_info
        if "MAIN-RELAX" in sb_file:
            try:
               with open('JARVIS-ID','r') as jf:
                     jl=jf.read()
                     name=jl.split('\n')[0] #str("JVASP-")+str(count)
            except:
               pass
            all_info['id']=name

            try:
               jfile=str(os.getcwd())+str("/")+str(sb_file)
               jdata = loadfn(jfile, cls=MontyDecoder)
               for data in jdata:
                   final_energy=float(data["final_energy"])/float(len(data["poscar_final"]))
                   final_str=(data["poscar_final"])
                   initial_pos=(data["poscar_initial"])
                   all_info['final_energy']=final_energy
                   all_info['initial_pos']=initial_pos
                   #all_info['composition']=final_str.composition
                   all_info['final_str']=final_str
            except:
               pass



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
              xphi=locpot(LOCPOTfile=LOCPOTfile,xml=v_xml,direction='x')
              yphi=locpot(LOCPOTfile=LOCPOTfile,xml=v_xml,direction='y')
              zphi=locpot(LOCPOTfile=LOCPOTfile,xml=v_xml,direction='z')
              phi=[xphi,yphi,zphi]
              #complete_dos = v_xml.complete_dos
              #mbandgap=float(complete_dos.get_gap())
            except:
               pass
            all_info['mbandgap']=mbandgap
            all_info['phi']=phi
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
            all_info['bbandgap']=bbandgap
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
            all_info['elastic']=elast
        if "Effective_Mass" in sb_file:
            try:
              efile=str(os.getcwd())+str("/")+str('Effective_Mass.json')
              eff= loadfn(efile, cls=MontyDecoder)
              n300=eff['n']['300'][0][0]
              p300=eff['p']['300'][0][0]
            except:
               pass
            ##all_info['eff']=eff
      os.chdir(pwd)
      c11=str(c11)+str(unit)
      all_info['dir']=dir
      #print "dir,phi=",dir,phi
      #print "bandgap=",mbandgap,bbandgap,c11,n300,p300
     
    except:
        pass
    return all_info    #mbandgap,bbandgap,c11,n300,p300,phi,opt_info

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
arr=[]
count=0
arg=0
for fil in glob.glob("*-*/*_*"):
#for fil in glob.glob("*PBEBO"):
    if 'PBEBO' in fil and '.png' not in fil:
        dir=str(os.getcwd())+str("/")+str(fil)
        print dir
        info=getprops(dir)
        arr.append(info)
fileXC=str('all_info3a')+'_data.json'
fjXC=open(fileXC,'w')

fjXC.write(json.dumps(arr,cls=MontyEncoder))
fjXC.close()

