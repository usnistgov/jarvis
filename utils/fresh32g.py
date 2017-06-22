#Note put printin ()
from monty.fractions import gcd
from monty.json import MontyEncoder
from collections import OrderedDict
from pymatgen.io.aseio import AseAtomsAdaptor
from pymatgen.io.vasp.outputs import Oszicar
from mpinterfaces.data_processor import MPINTVasprun
from pymatgen.matproj.rest import MPRester
import zipfile,sys
from pymatgen.io.vaspio_set import MPVaspInputSet, MPNonSCFVaspInputSet
from pymatgen.core.periodic_table import get_el_sp, Element
import requests,math,string
import fileinput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from bokeh.models import ColumnDataSource
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from ase.calculators.vasp import VaspChargeDensity
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn
from pymatgen.electronic_structure.core import Spin, Orbital
from bokeh.io import output_file, show, vform
from datetime import date
from random import randint
import copy,json
from bokeh.plotting import hplot
from bokeh.io import output_file, show, vform
from bokeh.models.widgets import TextInput
from bokeh.models import ColumnDataSource, OpenURL, TapTool
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import *
import periodic_table
#from bokeh.sampledata import periodic_table
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider,TapTool
from bokeh.plotting import figure, output_file, show
from bokeh.models.widgets import (
Select, DataTable, TableColumn, StringFormatter, HBox, VBox,
NumberFormatter, StringEditor, IntEditor, NumberEditor, SelectEditor)
import pandas as pd
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import vasp2boltz as vs
import os,shutil,time,datetime
from pymatgen.core.structure import Structure
import numpy as np
from asap3.analysis.rdf import RadialDistributionFunction
#from ase.lattice.compounds import *
from ase import data
import itertools
from pymatgen.electronic_structure  import boltztrap
from pymatgen.util.plotting_utils import get_publication_quality_plot
import matplotlib,yaml,os
from pymatgen.io.vasp.inputs import Potcar,Incar, Kpoints
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from pymatgen.electronic_structure.plotter import BSPlotterProjected, BSPlotter, DosPlotter
import matplotlib,time
matplotlib.use('Agg')
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import glob

MAPI_KEY = os.environ.get("MAPI_KEY", "")

def get_struct_from_mp(formula):
    ratio=0
    pf='na'
    compat=False
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        cnt=0

        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            mp=str(d['material_id'])
            bg=float(d["band_gap"])
            structure = m.get_structure_by_material_id(x['material_id'])
            a,b,c=structure.lattice.abc
            x['spacegroup'] = d['spacegroup']
            sg = x['spacegroup']['symbol']
            sg = sg.replace('/', 'slash')
            tm=str(d["total_magnetization"])
            pf=d["pretty_formula"]
            ucf=d["unit_cell_formula"]
            bg=d["band_gap"]
            ehull=d["e_above_hull"]
            icsd=d["icsd_id"]
            enp=d["energy_per_atom"]
#           print  x['spacegroup'],x['material_id']
            cnt=cnt+1
            filename=str(sg)+".POSCAR"
            line=str(filename)+'  '+'\n'
            structure = m.get_structure_by_material_id(x['material_id'],final=False)
            ini_structure = m.get_structure_by_material_id(x['material_id'],final=False)
            fin_structure = m.get_structure_by_material_id(x['material_id'])
            a1,b1,c1=structure.lattice.abc
            ratio= round(((float(c)-float(c1))/float(c1)),2)
            if ratio >=0.1 :
                print  x['material_id']
            if ratio >=0.1 and bg>0.1 and ehull==0 and icsd !=[] and icsd !=None:
               if  ehull==0 :
                print  "bandgap>0.1",x['material_id']
                compat=True
    return icsd #,ehull,tm,mp,pf,enp,bg,ini_structure,fin_structure
 
def ZipDir(inputDir, outputZip):
    '''Zip up a directory and preserve symlinks and empty directories'''
    zipOut = zipfile.ZipFile(outputZip, 'w', compression=zipfile.ZIP_DEFLATED)

    rootLen = len(os.path.dirname(inputDir))
    def _ArchiveDirectory(parentDirectory):
        #contents = os.listdir(parentDirectory)
        contents = ['INCAR', 'POSCAR', 'CONTCAR',  'KPOINTS']
        #store empty directories
        if not contents:
            #http://www.velocityreviews.com/forums/t318840-add-empty-directory-using-zipfile.html
            archiveRoot = parentDirectory[rootLen:].replace('\\', '/').lstrip('/')
            zipInfo = zipfile.ZipInfo(archiveRoot+'/')
            zipOut.writestr(zipInfo, '')
        for item in contents:
            fullPath = os.path.join(parentDirectory, item)
            if os.path.isdir(fullPath) and not os.path.islink(fullPath):
                _ArchiveDirectory(fullPath)
            else:
                archiveRoot = fullPath[rootLen:].replace('\\', '/').lstrip('/')
                if os.path.islink(fullPath):
                    # http://www.mail-archive.com/python-list@python.org/msg34223.html
                    zipInfo = zipfile.ZipInfo(archiveRoot)
                    zipInfo.create_system = 3
                    # long type of hex val of '0xA1ED0000L',
                    # say, symlink attr magic...
                    zipInfo.external_attr = 2716663808L
                    zipOut.writestr(zipInfo, os.readlink(fullPath))
                else:
                    zipOut.write(fullPath, archiveRoot, zipfile.ZIP_DEFLATED)
    _ArchiveDirectory(inputDir)

    zipOut.close()


def vtot(pref='',storedir=None,name=''):
#def vtot(lfile=None):
        for a in glob.glob(str('*')+str('.json')):
            if pref in a :#.startswith(pref):
		"""
		A script which averages a CHGCAR or LOCPOT file in one direction to make a 1D curve.
		User must specify filename and direction on command line.
		Depends on ase
		"""

                fold=a
                lfile=str(os.getcwd())+str("/")+str(a.split(".json")[0])+str("/")+str("LOCPOT")

		#if len(sys.argv) != 3:
		#    print "\n** ERROR: Must specify name of file and direction on command line."
		#    print "eg. vtotav.py LOCPOT z."                                            
		#    sys.exit(0)

		#if not os.path.isfile(sys.argv[1]):
		#    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
		#    sys.exit(0)

		# Read information from command line
		# First specify location of LOCPOT 
		LOCPOTfile = lfile#str(os.getcwd())+str("/")+str('LOCPOT') #sys.argv[1].lstrip()

		# Next the direction to make average in 
		# input should be x y z, or X Y Z. Default is Z.
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
                run=str(os.getcwd())+str("/")+str(fold.split(".json")[0])+str("/")+str("vasprun.xml")
               # run=str(os.getcwd())+str("/")+str(fold)+str("/vasprun.xml")
                vrun=Vasprun(run)
                Ef=vrun.efermi
                avg_max=max(average)
                dif=float(avg_max)-float(Ef)
                print ("Ef,max,wf=",Ef,avg_max,dif,os.getcwd())
		plt = get_publication_quality_plot(14, 10)
		plt.xlabel('z (Angstrom)')
		plt.plot(xs,ys,'-',linewidth=2,markersize=10)
                horiz_line_data = np.array([avg_max for i in xrange(len(xs))])
		plt.plot(xs,horiz_line_data,'-')
                horiz_line_data = np.array([Ef for i in xrange(len(xs))])
		plt.plot(xs,horiz_line_data,'-')
		plt.ylabel('Potential (eV)')
		ax = plt.gca()
		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		plt.title(str("Energy difference ")+str(round(float(dif),3))+str(" eV"),fontsize=26)
		filename=str(storedir)+str("/")+str(name)+str("/")+str('ESTAT.png')
		plt.tight_layout()
 
		plt.savefig(filename)
                plt.close()
		#print "\nEnd of calculation." 
                #print "Program was running for %.2f seconds." % runtime
                return dif





def periodic(html=None,csv=None,header=["Symbol","JV","REFID","En/at","Bv","Eg","eps0","gamma0","convg_kp","convg_encut"]):

		output_file(html)
		##FOR INFO TABLE
		#sym is content, row is period and column is group
		#p0 is the ELEMENTS ETC HEADER
		#p is the periodic table
		#p1 is the list
		#p2 is the warning
		TOOLS="tap,save"
		#TOOLS="resize,tap,save"
		#file='dat.json'

		sample=pd.read_csv(csv)



		period_range1=[str(x) for x in reversed(sorted(set(sample["period"])))]

		els=[]
		keys=[]
                max_period=1
		for i in range(0,len(sample)):
		    if sample["key"][i] not in keys:
                       group_count=1
                       period_count=1
		       keys.append(sample["key"][i])
		       els.append(sample["symbol"][i])
                    else:
                       group_count=group_count+1
                       period_count=period_count+1
                       if period_count>max_period:
                          max_period=period_count
                ##print group_count,max_period
                group_range1=[]
                group_count=len(header)
                for i in range(0,group_count):
                    group_range1.append(str(i+1))
                ##print group_range1,period_range1
		source_sample = ColumnDataSource(
		    data=dict(
			group=[str(x) for x in sample["group"]],
			period=[str(y) for y in sample["period"]],
			sym=sample["symbol"],
			symx=[str(x)+":0.1" for x in sample["group"]],
			els=els,
			keys=keys,
		    )
		)


                import sys
                #sys.exit()



		source_sample1 = ColumnDataSource(
		    data=dict(
			group=[],
			period=[],
			sym=[],
			symx=[],
			els=[],
			keys=[],
		    )
		)





		source_head = ColumnDataSource(
		    data=dict(
			group=[str(i+1) for i in range(1,len(header)+1)],
			period=[str(1) for i in range(1,len(header)+1)],
			#group=[str(i+1) for i in range(1,len(header)+1)],
			#period=[str(1) for i in range(1,len(header)+1)],
			sym=header,
		    )
		)
		per=[str(x) for x in reversed(sorted(set(source_head.data["period"])))]
		grp=[str(x) for x in range(2,len(header)+2)]
                print per
                print grp
                #sys.exit()
		TOOLS="resize,save"
		p0 = figure(responsive=True,
		    x_range=grp, y_range=per, tools=TOOLS)
		p0.rect("group", "period", 0.9, 0.9, source=source_head, fill_alpha=0.6)
		text_props = {
		    "source": source_head,
		    "angle": 0,
		    "text_align": "center",
		    "color": "black",
		    "text_baseline": "middle"
		}
 
                p0.border_fill_color = "whitesmoke"
		p0.min_border_top = 0
		p0.min_border_bottom = 0
		p0.plot_width = 900
		p0.plot_height = 90
		p0.logo = None
		p0.toolbar_location = None
		p0.grid.grid_line_color = None
		p0.axis.axis_line_color = None
		p0.axis.major_tick_line_color = None
		p0.axis.minor_tick_line_color = None
		p0.axis.major_tick_out = None
		p0.axis.major_tick_in = None
		p0.axis.major_label_text_font_size='0pt'  # turn off major ticks
		p0.axis.major_tick_line_color = None   # turn off minor ticks
		#p0.axis[1].ticker.num_minor_ticks = 0
		p0.text(x="group", y="period", text="sym",
		    text_font_style="bold", text_font_size="8pt", **text_props)

		print "grouprang1=",group_range1
		TOOLS="tap,save,resize"
		p1 = figure( tools=TOOLS,
		    x_range=group_range1, y_range=period_range1)
		p1.min_border_top = 0
                p1.border_fill_color = "whitesmoke"
		p1.min_border_bottom = 0
		p1.plot_width = 1200
		#p1.plot_width = 900
                print "max period=",max_period
		p1.plot_height = 600 #CHANGE HERE 
		#p1.toolbar_location = None
		p1.logo = None
		p1.toolbar_location = None
		p1.grid.grid_line_color = None
		p1.axis.axis_line_color = None
		p1.axis.major_tick_line_color = None
		p1.axis.minor_tick_line_color = None

		p1.axis.major_tick_out = None
		p1.axis.major_tick_in = None
		p1.toolbar_location = None
		p1.grid.grid_line_color = None
		p1.axis.axis_line_color = None
		p1.axis.major_label_text_font_size='0pt'  # turn off major ticks
		p1.axis.major_tick_line_color = None   # turn off minor ticks
		#p1.axis[1].ticker.num_minor_ticks = 0
		p1.rect("group", "period", 0.9, 0.9,
		    fill_alpha=0.6, color="type_color",source=source_sample1,name="foo")
		text_props = {
		    "source": source_sample1,
		    "angle": 0,
		    "color": "black",
		    "text_align": "left",
		    "text_baseline": "middle"
		}
		p1.text(x="symx", y="period", text="sym",
		    text_font_style="bold", text_font_size="8.0pt", **text_props)
		source_url = ColumnDataSource(
		    data=dict(
			url='',
		    )
		)
		#url = "https://www.materialsproject.org/materials/@sym/"
		taptool = p1.select(type=TapTool)
		renderer = p1.select(name="foo")[0]
		renderer.nonselection_glyph=renderer.glyph
		#renderer.nonselection_glyph=renderer.glyph.clone()
		#renderer.nonselection_glyph=renderer
		taptool.callback = CustomJS(args=dict(  source_sample1=source_sample1), code="""
			var inds = cb_obj.get('selected')['1d'].indices;
			var dat = source_sample1.get('data');
			inp=dat['sym'][inds[0]];
                        console.log(inds);
			console.log(("image.gif").startsWith("image"));
			if (inp.startsWith("mp-")){
		 
			url = "https://www.materialsproject.org/materials/"+inp
			}
			else if (inp.startsWith("JVASP-")){
			//url = "http://www.ctcms.nist.gov/~knc6/VASP_DFT/"+inp+".html"
			url = "http://www.ctcms.nist.gov/~knc6/jsmol/"+inp+".html"
			}

			else{
			url="www.ctcms.nist.gov/~knc6/periodic.html"
			}
			window.open(url, inp);
			//window.open(url, 'blank');
			//window.open(url, '_blank');
			//return false
		    """)
		#taptool.callback=OpenURL(url="@url")
		#taptool = p1.select(type=TapTool, action=OpenURL(url=source_url.data['url']))

		#url="@url"
		#taptool.callback = OpenURL(url=url)



		callback_sample = CustomJS(args=dict( source_sample=source_sample, source_sample1=source_sample1), code="""
			var f = cb_obj.get('value') 
			//fin=[str("-".join(e)) for e in list(itertools.permutations(el.split("-")))]
			console.log(f);
			var f=f.split("-").sort();
			var f=f.join("-");
			console.log(f);
			var data = source_sample.get('data');
			var data1 = source_sample1.get('data');
			sym=data["sym"];
			symx=data["symx"];
			grp=data["group"];
			per=data["period"];
			max=grp[0];
			for (i=0;i<grp.length;i++){
			    if(grp[i]>max){
			       max=grp[i];
			       
			 }
			}
			data1["group"]=[];
			data1["period"]=[];
			data1["sym"]=[];
			data1["symx"]=[];
			for (i=0;i<sym.length;i++){
			   if (sym[i]==f){
			     for (j=0;j<11;j++){
			  // for (j=0;j<16;j++){  CHANGE NUMBER HERE
			       data1['group'].push(data['group'][i+j]);
			       data1['period'].push(data['period'][i+j]);
			       data1['sym'].push(data['sym'][i+j]);
			       data1['symx'].push(data['symx'][i+j]);
			   
			 }
			     }
			}

			source_sample1.trigger('change');
		    """)

		text_input1 = TextInput(value="Enter", title=str("SEARCH : "),callback=callback_sample)
		#text_input1 = TextInput(value="Search: ", title=str("Available : ")+str(available),callback=callback_sample)

		###	FOR PERIODIC TABLE ##
		periodic_table.elements["atomic mass"] = periodic_table.elements["atomic mass"].astype(str)

		elements = periodic_table.elements[periodic_table.elements["group"] != "-"]

		group_range = [str(x) for x in range(1,19)]
		period_range = [str(x) for x in reversed(sorted(set(elements["period"])))]

		periodic_table.elements["atomic mass"] = periodic_table.elements["atomic mass"].astype(str)

		elements = periodic_table.elements[periodic_table.elements["group"] != "-"]

		group_range = [str(x) for x in range(1,19)]
		period_range = [str(x) for x in reversed(sorted(set(elements["period"])))]

		colormap = {
		    "alkali metal"         : "#a6cee3",
		    "alkaline earth metal" : "#1f78b4",
		    "halogen"              : "#fdbf6f",
		    "metal"                : "#b2df8a",
		    "metalloid"            : "#33a02c",
		    "noble gas"            : "#bbbb88",
		    "nonmetal"             : "#baa2a6",
		    "transition metal"     : "#e08e79",
		    "lanthanoid"     : "#a39d33",
		    "actinoid"     : "#a37333",
		}

		source = ColumnDataSource(
		    data=dict(
			group=[str(x) for x in elements["group"]],
			period=[str(y) for y in elements["period"]],
			symx=[str(x)+":0.1" for x in elements["group"]],
			numbery=[str(x)+":0.8" for x in elements["period"]],
			massy=[str(x)+":0.15" for x in elements["period"]],
			namey=[str(x)+":0.3" for x in elements["period"]],
			sym=elements["symbol"],
			name=elements["name"],
			cpk=elements["CPK"],
			atomic_number=elements["atomic number"],
			electronic=elements["electronic configuration"],
			mass=elements["atomic mass"],
			type=elements["metal"],
			type_color=[colormap[x] for x in elements["metal"]],
		    )
		)


		TOOLS = "tap,save,hover,resize"
		#TOOLS = "resize,tap,save"

		p = figure(title="JARVIS for DFT (2D Materials)", tools=TOOLS,
		#p = figure( tools=TOOLS,
		    x_range=group_range, y_range=period_range)
		p.plot_width = 1200
		p.logo = None
		p.toolbar_location = None
		#p.toolbar_location = "left"
		p.min_border_top = 0
		p.min_border_bottom = 0
		#p.plot_width = 1200
		p.plot_height = 600
		p.logo = None
		p.toolbar_location = None
		p.grid.grid_line_color = None
		p.axis.axis_line_color = None
		p.axis.major_tick_line_color = None
		p.axis.minor_tick_line_color = None
		p.axis.major_tick_out = None
		p.axis.major_tick_in = None
		p.axis.major_label_text_font_size='0pt'  # turn off major ticks
		p.axis.major_tick_line_color = None   # turn off minor ticks
		#p.axis[1].ticker.num_minor_ticks = 0

		p.rect("group", "period", 0.9, 0.9, source=source,
		    fill_alpha=0.6, color="type_color",name="foo1")

		text_props = {
		    "source": source,
		    "angle": 0,
		    "color": "black",
		    "text_align": "left",
		    "text_baseline": "middle"
		}

		p.text(x="symx", y="period", text="sym",
		    text_font_style="bold", text_font_size="13pt", **text_props)

		p.text(x="symx", y="numbery", text="atomic_number",
		    text_font_size="7pt", **text_props)

		p.text(x="symx", y="namey", text="name",
		    text_font_size="7pt", **text_props)

		p.text(x="symx", y="massy", text="mass",
		    text_font_size="5pt", **text_props)

		p.grid.grid_line_color = None

		hover = p.select(dict(type=HoverTool))
		hover.tooltips = OrderedDict([
		    ("name", "@name"),
		    ("atomic number", "@atomic_number"),
		    ("type", "@type"),
		#    ("atomic mass", "@mass"),
		#    ("CPK color", "$color[hex, swatch]:cpk"),
		    ("electronic configuration", "@electronic"),
		])



		#plot = figure(plot_width=400, plot_height=400)


		tap = p.select(dict(type=TapTool))
		#tap.callback = CustomJS(code='''alert("pressed");''')
		tap.callback = CustomJS(args=dict( source=source, source_sample=source_sample,source_sample1=source_sample1), code="""
			var inds = cb_obj.get('selected')['1d'].indices;
			var dat = source.get('data');
			inp=[dat['sym'][inds[0]]];
			for (i=1 ; i<inds.length; i++){ 
			inp=inp+"-"+(dat['sym'][inds[i]])
			console.log(inp);
			var inp=inp.split("-").sort()
			var inp=inp.join("-")
			//console.log(inp);
			}
			//console.log(inds,inp)
			//fin=[str("-".join(e)) for e in list(itertools.permutations(el.split("-")))]
			var data = source_sample.get('data');
			var data1 = source_sample1.get('data');
			sym=data["sym"];
			symx=data["symx"];
			per=data["period"];
			grp=data["group"];
			max=grp[0];
			for (i=0;i<grp.length;i++){
			    if(grp[i]>max){
			       max=grp[i];
			 }
			}
			//console.log(max)
			data1["group"]=[];
			data1["period"]=[];
			data1["sym"]=[];
			data1["symx"]=[];
			for (i=0;i<sym.length;i++){
			   if (sym[i]==inp){
			     for (j=0;j<11;j++){
			  // for (j=0;j<16;j++){ #CHANGE NUMBER HERE
			       data1['group'].push(data['group'][i+j]);
			       data1['period'].push(data['period'][i+j]);
			       data1['sym'].push(data['sym'][i+j]);
			       data1['symx'].push(data['symx'][i+j]);
			   
			 }
			     }
			}

			source_sample1.trigger('change');
		    """)

		renderer = p.select(name="foo1")[0]
		###renderer.nonselection_glyph=renderer.glyph.clone()



		layout = vform(text_input1,p,p0, p1)
		#layout = vform(text_input1,p,p2,p0, p1)
		save(layout)
		#show(layout)
                f=open(html,'a')
                add=str('<body style="background-color:lightgreen;">')


                for line in fileinput.FileInput(html,inplace=1):
                    if '<html lang="en">' in line:
                        line=line.replace(line,line+add)
                    print line,
                f.close()
		#show(p)
#periodic(html='ff.html',header=['CALCS','NN','MM','BB','MPID','FORMULA','EN/ATOM (eV)','C11(GPa)','C22(GPa)','C33(GPa)','C12(GPa)','C13(GPa)','C23(GPa)','C44(GPa)','Bv(GPa)','Gv(GPa)','E_hull_mp','FORCEFIELD'],csv='EAM_alloy16.csv')
plt = get_publication_quality_plot(14, 10)
plt.tight_layout()
def boltz_bader(pref='',storedir=None,name=''):
        from pymatgen.io.vaspio.vasp_output import Vasprun
        from pymatgen.electronic_structure.plotter import BSPlotter
        import vasp2boltz as vs
        from pymatgen.electronic_structure  import boltztrap
        from pymatgen.electronic_structure.boltztrap  import BoltztrapAnalyzer,BoltztrapPlotter
        print ("CWD=",os.getcwd())
        for a in glob.glob(str('*')+str('.json')):
            print ("a,pref",a,pref)
            if pref in a :#.startswith(pref):
            #if a.startswith(pref):
                print ("a,pref boltz",a,pref)
                fold=a
                run=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                vrun=Vasprun(run)
                print ("vrun=",vrun)
                fin_en= vrun.final_energy
                kpfile=str(a.split(".json")[0])+str("/")+str("KPOINTS")
                outfile=str(a.split(".json")[0])+str("/")+str("OUTCAR")
                contcar=Structure.from_file(str(a.split(".json")[0])+str("/")+str("CONTCAR"))
                print ("FIN ENNNN=",fin_en,type(fin_en),contcar.composition.num_atoms)
                target_file=str(str(a.split(".json")[0]))
                dest=str(storedir)+str("/")+str(name)+str("/")+str('MAIN-RELAX')+str(".zip")
                ZipDir(target_file,dest)
                outf_file=open(outfile,'r')
                lines = outf_file.read().splitlines()
                PSCENC='na'
                TEWEN='na'
                DENC='na'
                EXHF='na'
                XCENC='na'
                PAW_DB='na'
                EENTRO='na'
                EATOM='na'
                TOTEN='na'
                additional_info=[]
                for i,l in enumerate(lines):
                       
                    if "Free energy of the ion-electron system (eV)" in l:
                         additional_info=[]
                         #PSCENC=str(lines[i+2].split('=')[1]).split()[0]#lines[i+2]
                         PSCENC=str(lines[i+2].split('=')[1])#lines[i+2]
                         additional_info.append(PSCENC)       
                         TEWEN=str(lines[i+3].split('=')[1])#.split()[0]#lines[i+3]
                         additional_info.append(TEWEN)       
                         DENC=str(lines[i+4].split('=')[1])#.split()[0]#lines[i+4]
                         additional_info.append(DENC)       
                         EXHF=str(lines[i+5].split('=')[1])#.split()[0]#lines[i+5]
                         additional_info.append(EXHF)       
                         XCENC=str(lines[i+6].split('=')[1])#.split()[0]#lines[i+6]
                         additional_info.append(XCENC)       
                         PAW_DB=str(lines[i+7].split('=')[1])#.split()[0]#lines[i+7]
                         additional_info.append(PAW_DB)       
                         EENTRO=str(lines[i+8].split('=')[1])#.split()[0]#lines[i+8]
                         additional_info.append(EENTRO)       
                         EBANDS=str(lines[i+9].split('=')[1])#.split()[0]#lines[i+9]
                         additional_info.append(EBANDS)       
                         EATOM=str(lines[i+10].split('=')[1])#.split()[0]#lines[i+10]
                         additional_info.append(EATOM)       
                         #TOTEN=str(lines[i+13]).split('=')[1]#.split('=')[1])#.split()[0]).split('eV')[0]#lines[i+12]
                         #additional_info.append(TOTEN)       
                outf_file.close()
                outf_file=open(outfile,'r')
                found_mag=False
                mag_mom='na'
                nelect='na'
                for i,l in enumerate(outf_file):
                    if "NELECT" in l:
                       nelect=str(l.split("=")[1])
                    if "magnetization (x)" in l:
                         found_mag=True
                    if 'tot' in l and "#" not in l and found_mag==True and 'total' not in l:
                          mag_mom=str(l.split()[4])
                          print ("MAGNETIC MOMENT=",mag_mom)
                          found_mag=False
                additional_info.append(mag_mom)       
                additional_info.append(nelect)       
                outf_file.close()
                #additional_info.append([PSCENC,TEWEN, DENC,EXHF,XCENC,PAW_DB,EENTRO,EATOM,mag_mom])
                print ("additional info=", additional_info )
        #print ("DFGHJKL;DFGHJKL",a)
        #time.sleep(5)
        #fin_en=float(fin_en)/float(contcar.composition.num_atoms)
        el_list=sorted(list(contcar.symbol_set))      
        search= ('-'.join([item for item in el_list]))  
        finder = SpacegroupAnalyzer(contcar)
        ##num=finder.get_spacegroup_number()
        num=finder.get_spacegroup_symbol()
        num.replace('/', 'slash')
        mat_cvn = finder.get_conventional_standard_structure()
#      PLOT XRD
        c = XRDCalculator()
        filename=str(storedir)+str("/")+str(name)+str("/")+str('xrd.png')
        c.get_xrd_plot(contcar,annotate_peaks=False).savefig(filename)
        #c.get_xrd_plot(contcar,annotate_peaks=False).tight_layout().savefig(filename)
        c.get_xrd_plot(contcar,annotate_peaks=False).close()
#      PLOT RDF
        #try:
	rng=40.0
	bins = 200
        print ("NOW CONTCAR FOR ATOMS IS",contcar)
	atoms= AseAtomsAdaptor().get_atoms(contcar)
        try:		
           sa1=int(float(rng)/float(max(abs(atoms.get_cell()[0]))))+2
           sa2=int(float(rng)/float(max(abs(atoms.get_cell()[1]))))+2
           sa3=int(float(rng)/float(max(abs(atoms.get_cell()[2]))))+2
           print ('scal=',sa1,sa2,sa3)
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
        rng=15.0
	x = np.arange(bins) * rng / bins
	RDFobj = RadialDistributionFunction(atoms, rng, bins, verbose=True)
	symb_comb=[]
        plt = get_publication_quality_plot(14, 10)
	print ("symbs",symbs)
	#plt.close()
	if len(symbs)>1:
	   for el1 in  symbs:
	      for el2 in  symbs:
		 el=str(el1)+str('-')+str(el2)
		 if str(el1)+str('-')+str(el2) not in symb_comb and str(el2)+str('-')+str(el1) not in symb_comb:
		     symb_comb.append(el)
		     rdf = RDFobj.get_rdf(elements=(data.atomic_numbers[str(el1)], data.atomic_numbers[str(el2)]))
                     plt.tight_layout()
		     plt.plot(x, rdf,label=str(el1)+str('-')+str(el2),linewidth=2)
	elif len(symbs)==1:
		     rdf = RDFobj.get_rdf(elements=(data.atomic_numbers[str(symbs[0])], data.atomic_numbers[str(symbs[0])]))
        #             plt = get_publication_quality_plot(14, 10)
                     plt.tight_layout()
		     plt.plot(x, rdf,label=str(symbs[0])+str('-')+str(symbs[0]),linewidth=2)

        #plt = get_publication_quality_plot(14, 10)
	plt.legend(prop={'size':26})
	plt.xlabel("r (Angstrom)")
	plt.ylabel("g(r)")
        plt.tight_layout()
	filename=str(storedir)+str("/")+str(name)+str("/")+str('rdf.png')
	plt.savefig(filename)
	plt.close()
        #except:
        #      pass


        c_size=12
        dim1=int((float(c_size)/float( max(abs(contcar.lattice.matrix[0])))))+1
        dim2=int(float(c_size)/float( max(abs(contcar.lattice.matrix[1]))))+1
        dim3=int(float(c_size)/float( max(abs(contcar.lattice.matrix[2]))))+1
        super_contcar=contcar
        super_contcar.make_supercell([dim1,dim2,dim3])
        tmp1=AseAtomsAdaptor().get_atoms(super_contcar)
        try:
           tmp1.center(vacuum=0.0, axis=2)
           print ("tmp111111111111",tmp1)
           if tmp1.get_cell()[2][2]==0.0:
              print "in get cell"
              tmp1.set_cell([(tmp1.get_cell()[0]),(tmp1.get_cell()[1]),(tmp1.get_cell()[2][0],tmp1.get_cell()[2][1],3.0)])
              print ("tmpxxxxx",tmp1)
        except:
           #tmp1.center(vacuum=1.2, axis=2)
           pass
        tmp2=AseAtomsAdaptor().get_structure(tmp1)
        filename=str(storedir)+str("/")+str(name)+str("/")+str(name)+str('_prem.cif')
        tmp2.to(fmt= "cif", filename= filename)
        #contcar.to(fmt= "cif", filename= filename)
        #super_contcar.to(fmt= "cif", filename= filename)
        #contcar.to(fmt= "cif", filename= filename)
        #mat_cvn.to(fmt= "cif", filename= filename)
        avg=[[]]
	v = Vasprun(run, parse_projected_eigen = True)
	path=str(os.getcwd())+str("/")+str(fold.split(".json")[0])
	#path=str(a.split(".json")[0])
	print ("path a",path,a)
	v = Vasprun(run)
	#v = Vasprun(run, parse_projected_eigen = True)
	bs = v.get_band_structure(kpfile, line_mode = False)
	for line in open(outfile, 'r'):
	     if 'NELECT' in line:
		 nelec = float(line.split()[2])
	print (nelec)
        locf=str(path)+str("/LOCPOT")
        #vtot(lfile=locf)
	folder=str(path)+str('/')+str('boltztrap')
        #st = os.stat(folder)
        #skip=False
	#cs=boltztrap.BoltztrapRunner(bs,nelec)
        #try:
        #   os.remove(folder)
        #except:
        #  pass
	#cs.run(path_dir=path)
        if not os.path.exists(folder):
	   cs=boltztrap.BoltztrapRunner(bs,nelec)
	   cs.run(path_dir=path)
           #while skip!=True:
           #    if time.time() - st.st_mtime > timeout:
           #       errors.append('Frozenjob')





        try:

			v = Vasprun(run, parse_projected_eigen = True)
			path=str(os.getcwd())+str("/")+str(fold.split(".json")[0])
			#path=str(a.split(".json")[0])
			print ("path a",path,a)
			v = Vasprun(run)
			#v = Vasprun(run, parse_projected_eigen = True)
			bs = v.get_band_structure(kpfile, line_mode = False)
			for line in open(outfile, 'r'):
			     if 'NELECT' in line:
				 nelec = float(line.split()[2])
			print (nelec)
			folder=str(path)+str('/')+str('boltztrap')
			#st = os.stat(folder)
			#skip=False
			cs=boltztrap.BoltztrapRunner(bs,nelec)
			#try:
			#   os.remove(folder)
			#except:
			#  pass
			#cs.run(path_dir=path)
			if not os.path.exists(folder):
			   cs=boltztrap.BoltztrapRunner(bs,nelec)
			   cs.run(path_dir=path)
			   #while skip!=True:
			   #    if time.time() - st.st_mtime > timeout:
			   #       errors.append('Frozenjob')

			print ("SEE RUNNING IN",folder,'  ',name)
			bs = BoltztrapAnalyzer.from_files(folder)
			print ("SEE RUNNING IN2",folder)
			f = open('Effective_Mass','w')
			eff_m = open('Effective_Mass.json','w')
			eff_m.write(json.dumps(bs.get_average_eff_mass(),indent=4))
			eff_m.close() #write(json.dumps(bs.get_average_eff_mass(),indent=4)
			avg=bs.get_average_eff_mass(output='tensor')
			#line=str(bs.get_average_eff_mass()['p'][300])+'\n'
			line=str(bs.get_average_eff_mass())+'\n'
			f.write(line)
			f.close()
			plt = get_publication_quality_plot(14, 10)
			bp = BoltztrapPlotter(bs)
			cond_mu1=bp.plot_conductivity_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('Cond_mu.png')
			print ("cond filename",filename)
			cond_mu1.tight_layout()
			cond_mu1.savefig(filename)
			cond_mu1.close()
			#cond_mu.tight_layout().savefig(filename)
                        try:
                          plt.close()
                        except:
                             pass
			plt = get_publication_quality_plot(14, 10)
			zt_mu1=bp.plot_zt_mu(temp=100)
			#zt_mu=bp.plot_zt_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('ZT_mu_100.png')
			#plt.tight_layout()
			zt_mu1.tight_layout()
			zt_mu1.savefig(filename)
			zt_mu1.close()
                        try:
                          plt.close()
                        except:
                             pass
			#zt_mu.tight_layout().savefig(filename)
		       


			plt = get_publication_quality_plot(14, 10)
			zt_mu2=bp.plot_zt_mu(temp=300)
			#zt_mu=bp.plot_zt_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('ZT_mu_300.png')
			zt_mu2.tight_layout()
			zt_mu2.savefig(filename)
			zt_mu2.close()
			#zt_mu.tight_layout().savefig(filename)
			#plt.close()
                        try:
                          plt.close()
                        except:
                             pass


			plt = get_publication_quality_plot(14, 10)
			zt_mu3=bp.plot_zt_mu(temp=300,relaxation_time=1e-12)
			#zt_mu=bp.plot_zt_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('ZT_mu_300_rel_12.png')
			zt_mu3.tight_layout()
			zt_mu3.savefig(filename)
			zt_mu3.close()
                        try:
                          plt.close()
                        except:
                             pass
			#zt_mu.tight_layout().savefig(filename)
			#plt.close()



			plt = get_publication_quality_plot(14, 10)
			zt_mu4=bp.plot_zt_mu(temp=300,relaxation_time=1e-16)
			#zt_mu=bp.plot_zt_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('ZT_mu_300_rel_16.png')
			zt_mu4.tight_layout()
			zt_mu4.savefig(filename)
			zt_mu4.close()
			#zt_mu.tight_layout().savefig(filename)
			#plt.close()

                        try:
                          plt.close()
                        except:
                             pass



			plt = get_publication_quality_plot(14, 10)
			zt_mu5=bp.plot_zt_mu(temp=600)
			#zt_mu=bp.plot_zt_mu()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('ZT_mu_600.png')
			zt_mu5.tight_layout()
			zt_mu5.savefig(filename)
			#zt_mu.tight_layout().savefig(filename)
			#plt.close()
			zt_mu5.close()
                        try:
                          plt.close()
                        except:
                             pass

			plt = get_publication_quality_plot(14, 10)
			s_mu1=bp.plot_seebeck_mu(temp=100)
			s_mu1.tight_layout()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('S_mu_100.png')
			s_mu1.savefig(filename)
			s_mu1.close()

                        try:
                          plt.close()
                        except:
                             pass

			plt = get_publication_quality_plot(14, 10)
			s_mu2=bp.plot_seebeck_mu(temp=300)
			s_mu2.tight_layout()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('S_mu_300.png')
			s_mu2.savefig(filename)
			s_mu2.close()
                        try:
                          plt.close()
                        except:
                             pass


			plt = get_publication_quality_plot(14, 10)
			s_mu3=bp.plot_seebeck_mu(temp=600)
			s_mu3.tight_layout()
			filename=str(storedir)+str("/")+str(name)+str("/")+str('S_mu_600.png')
			s_mu3.savefig(filename)
	                s_mu3.close()
                        try:
                          plt.close()
                        except:
                             pass
        except:
          pass

        #fin_en=round(float(fin_en),4)
        return fin_en,search,num,avg,additional_info
def K_fig(pref='',storedir=None):
    x=[]
    y=[]
    for a in glob.glob(str('*')+str('.json')):
            if a.startswith('MAIN-RELAX'):
               main_kp=str(a.split(".json")[0])+str("/")+str("KPOINTS")
               main_kp_obj=Kpoints.from_file(main_kp)
               [k1,k2,k3]=main_kp_obj.kpts[0]
              # print ("kpoints====",k1,k2,k3)
            elif a.startswith(pref):
                k=int(float(str(a.split('-')[-1]).split('.json')[0]))
                #k=int(float(str(a.split('-')[1]).split('.json')[0]))
                #run=str(a.split(".json")[0])+str("/")+str("OSZICAR")
                run=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                vrun=Vasprun(run)
                #vrun=Oszicar(run)

                kpfile=str(a.split(".json")[0])+str("/")+str("KPOINTS")
                contcar=Structure.from_file((a.split(".json")[0])+str("/")+str("CONTCAR"))
                kpoints=Kpoints.from_file(kpfile)
                [xx,yy,zz]=kpoints.kpts[0]
                en =float(vrun.final_energy)#/float(contcar.composition.num_atoms)
                #en =float(vrun.final_energy)/float(contcar.composition.num_atoms)
               # print "ENERGYYY,AT",en,float(contcar.composition.num_atoms)
                x.append(k)
                y.append(en)
    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    len_xs=len(xs)
    xs1=[]
    ys1=[]
    target_ys=ys[-1]
    #print "target=",target_ys
    for i,el in enumerate(ys):
        if el <= (float(target_ys)+0.002) :
            xs1.append(xs[i])
            ys1.append(ys[i])
    #print "xs,ys=",xs,ys
    #print "xs1,ys1=",xs1,ys1
    left, bottom, width, height = [0.5, 0.5, 0.35, 0.35]
    plt = get_publication_quality_plot(14, 10)
    fig, ax1 = plt.subplots()
    plt.xlabel('Increment in K point',fontsize=20)
    plt.ylabel('Energy (eV)',fontsize=20)
    plt.title(str("Converged at ")+str(k1)+str("x")+str(k2)+str("x")+str(k3)+str(" ")+str("Automatic Mesh ")+str(int(k)-25),fontsize=26)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax2 = fig.add_axes([left, bottom, width, height])
    #ax2.plot(xs1,ys1, '.-',linewidth=2,markersize=10)
    ax1.plot(xs,ys, 's-',linewidth=2,markersize=10)
    el_list=sorted(list(contcar.symbol_set))
    search= ('-'.join([item for item in el_list]))
    #ax1.xlabel('Increment in K point')
    #print "xs,ys"
    #print xs
    #print ys
    plt.plot(xs1,ys1,'.-',linewidth=2,markersize=10)
    plt.ylim([float(target_ys)+0.002,float(target_ys)-0.002])
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    filename=str(storedir)+str("/")+str('KDen.png')
    k_convg=str(xx)+str("x")+str(yy)+str("x")+str(zz)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    #print("search2",search,k_convg)
    time.sleep(5)
    return  search,k_convg

def K_figold(pref='',storedir=None):
    x=[]
    y=[]
    for a in glob.glob(str('*')+str('.json')):
            if a.startswith('MAIN-RELAX'):
               main_kp=str(a.split(".json")[0])+str("/")+str("KPOINTS") 
               main_kp_obj=Kpoints.from_file(main_kp)
               [k1,k2,k3]=main_kp_obj.kpts[0]
               print "kpoints====",k1,k2,k3
            elif a.startswith(pref):
                k=int(float(str(a.split('-')[-1]).split('.json')[0]))
                #k=int(float(str(a.split('-')[1]).split('.json')[0]))
                #run=str(a.split(".json")[0])+str("/")+str("OSZICAR")
                run=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                vrun=Vasprun(run)
                #vrun=Oszicar(run)
          
                kpfile=str(a.split(".json")[0])+str("/")+str("KPOINTS")
                contcar=Structure.from_file((a.split(".json")[0])+str("/")+str("CONTCAR"))
                kpoints=Kpoints.from_file(kpfile)
                [xx,yy,zz]=kpoints.kpts[0]
                en =float(vrun.final_energy)#/float(contcar.composition.num_atoms)
                #en =float(vrun.final_energy)/float(contcar.composition.num_atoms)
                print "ENERGYYY,AT",en,float(contcar.composition.num_atoms)
                x.append(k)
                y.append(en)
    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    el_list=sorted(list(contcar.symbol_set))      
    search= ('-'.join([item for item in el_list]))  
    plt = get_publication_quality_plot(14, 10)
    plt.xlabel('Increment in K point')
    plt.ylabel('Energy (eV)')
    print "xs,ys"
    print xs
    print ys
    plt.plot(xs,ys,'s-',linewidth=2,markersize=10)
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.title(str("Converged at ")+str(k1)+str("x")+str(k2)+str("x")+str(k3)+str(" ")+str("Automatic Mesh ")+str(int(k)-25),fontsize=26)
    filename=str(storedir)+str("/")+str('KDen.png')
    k_convg=str(xx)+str("x")+str(yy)+str("x")+str(zz)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print("search2",search,k_convg)
    time.sleep(5)
    return  search,k_convg
def EN_fig(pref='',storedir=None):
    x=[]
    y=[]
    for a in glob.glob(str('*')+str('.json')):
            if a.startswith('MAIN-RELAX'):
               
               main_inc=str(a.split(".json")[0])+str("/")+str("INCAR") 
               main_inc_obj=Incar.from_file(main_inc)
               
               convg_encut=float(main_inc_obj['ENCUT'])
               #print ("IN EN_fig",main_inc)
            elif a.startswith(pref):
                e=int(str(a.split('-')[-1]).split('.json')[0])
                #e=int(str(a.split('-')[1]).split('.json')[0])
                run=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                contcar=Structure.from_file((a.split(".json")[0])+str("/")+str("CONTCAR"))
                #vrun=MPINTVasprun(run)
                #run=str(a.split(".json")[0])+str("/")+str("OSZICAR")
                #vrun=Oszicar(run)#Vasprun(run, parse_dos=False, parse_potcar_file=False)
                vrun=Vasprun(run)
                infile=str(a.split(".json")[0])+str("/")+str("INCAR")
                inc=Incar.from_file(infile)
                encut=inc['ENCUT']
                en =float(vrun.final_energy)#/float(contcar.composition.num_atoms)
                #en =float(vrun.final_energy)/float(contcar.composition.num_atoms)
                x.append(encut)
                y.append(en)
    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    plt = get_publication_quality_plot(14, 10)
    plt.ylabel('Energy (eV)')
    plt.plot(xs,ys,'s-',linewidth=2,markersize=10)
    plt.xlabel('Increment in ENCUT ')
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.title(str("Converged at ")+str(int(convg_encut))+str("eV"),fontsize=26)
    filename=str(storedir)+str("/")+str('Encut.png')
    plt.tight_layout()
    plt.savefig(filename)
    enc_convg=str(encut)
    #plt.savefig(str(pref)+str(".png"))
    plt.close()
def optics(pref='',storedir=None):
    for a in glob.glob(str('*')+str('.json')):
            if a.startswith(pref):
                ru=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                try:
                        contcar=Structure.from_file(str(a.split(".json")[0])+str("/")+str("POSCAR"))
                        if "Surf" in a:
                            print "yes optic",a
                            ratio_c=0.1*float(abs(contcar.lattice.matrix[2][2]))#*(10**9)*(10**-10) #N/m unit
                except:
                        pass
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





 
                #plt.close()
                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,realx,linewidth=2,label=r'$\epsilon_1x$')
                plt.plot(en,realy,linewidth=2,label=r'$\epsilon_1y$')
                plt.plot(en,realz,linewidth=2,label=r'$\epsilon_1z$')
                plt.legend(prop={'size':26})
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Dielec. function (Real part)')
                filename=str(storedir)+str("/")+str('Real.png')
                plt.tight_layout()
                plt.savefig(filename)
                plt.tight_layout()
                plt.close()
                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,imagx,linewidth=2,label=r'$\epsilon_2x$')
                plt.plot(en,imagy,linewidth=2,label=r'$\epsilon_2y$')
                plt.plot(en,imagz,linewidth=2,label=r'$\epsilon_2z$')
                plt.xlabel('Energy (eV)')
                plt.xlim([0,30])
                plt.ylabel('Dielec. function(Imag part)')
                filename=str(storedir)+str("/")+str('Imag.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.close()



                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,refrx,linewidth=2,label=r'$n_x$')
                plt.plot(en,refry,linewidth=2,label=r'$n_y$')
                plt.plot(en,refrz,linewidth=2,label=r'$n_z$')
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Refractive index')
                filename=str(storedir)+str("/")+str('Refr.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.close()


                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,extcx,linewidth=2,label=r'$k_x$')
                plt.plot(en,extcy,linewidth=2,label=r'$k_y$')
                plt.plot(en,extcz,linewidth=2,label=r'$k_z$')
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Extinction Coefficient')
                filename=str(storedir)+str("/")+str('Extc.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.tight_layout()
                plt.close()



                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,absorpx,linewidth=2,label=r'$\alpha_x$')
                plt.plot(en,absorpy,linewidth=2,label=r'$\alpha_y$')
                plt.plot(en,absorpz,linewidth=2,label=r'$\alpha_z$')
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Absorption coefficient')
                filename=str(storedir)+str("/")+str('Absorp.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.tight_layout()
                plt.close()


                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,eelsx,linewidth=2,label=r'$e_x$')
                plt.plot(en,eelsy,linewidth=2,label=r'$e_y$')
                plt.plot(en,eelsz,linewidth=2,label=r'$e_z$')
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Energy loss spectrum')
                filename=str(storedir)+str("/")+str('ELS.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.tight_layout()
                plt.close()



                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,opt_conx,linewidth=2,label=r'$\sigma_x$')
                plt.plot(en,opt_cony,linewidth=2,label=r'$\sigma_y$')
                plt.plot(en,opt_conz,linewidth=2,label=r'$\sigma_z$')
                plt.xlim([0,30])
                plt.xlabel('Energy (eV)')
                plt.ylabel('Optical conductivity (Real part)')
                filename=str(storedir)+str("/")+str('Opt_con.png')
                plt.tight_layout()
                plt.legend(prop={'size':26})
                plt.savefig(filename)
                plt.tight_layout()
                plt.close()


def banddosHSE(pref='',storedir=None):
    dirs_list=[d for d in os.listdir('./') if os.path.isdir(os.path.join('./', d))]
    print "sirlist=",dirs_list
    for l in dirs_list:
            print l
    #for a in glob.glob(str('*')+str('.json')):
            if l.startswith(pref):
               # from pymatgen.io.vasp.outputs import Vasprun
                from pymatgen.io.vaspio.vasp_output import Vasprun
                ru=str(l)+str("/")+str("vasprun.xml")
                kpfile=str(l)+str("/")+str("KPOINTS")




                run = Vasprun(ru, parse_projected_eigen = True)
                bands = run.get_band_structure(kpfile, line_mode = True, efermi = run.efermi)
                bsp =  BSPlotter(bands)
                zero_to_efermi=True
                data=bsp.bs_plot_data(zero_to_efermi)
                from pymatgen.util.plotting_utils import get_publication_quality_plot
                plt = get_publication_quality_plot(12, 8)
    #plt = get_publication_quality_plot(12, 8)
                band_linewidth = 3
                x_max = data['distances'][-1][-1]
                #print (x_max)
                for d in range(len(data['distances'])):
                   for i in range(bsp._nb_bands):
                      plt.plot(data['distances'][d],
                             [data['energy'][d]['1'][i][j]
                              for j in range(len(data['distances'][d]))], 'b-',
                             linewidth=band_linewidth)
                      if bsp._bs.is_spin_polarized:
                         plt.plot(data['distances'][d],
                                 [data['energy'][d]['-1'][i][j]
                                  for j in range(len(data['distances'][d]))],
                                 'r--', linewidth=band_linewidth)
                bsp._maketicks(plt)
                if bsp._bs.is_metal():
                     e_min = -10
                     e_max = 10
                     band_linewidth = 3

                for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)

                        for vbm in data['vbm']:
                            plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                        s=100)


        # Main X and Y Labels
                plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
                ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi \
                   else r'$\mathrm{Energy\ (eV)}$'
                plt.ylabel(ylabel, fontsize=30)
                plt.ylim(-5,10)
#    plt.ylim(data['vbm'][0][1] + e_min,
#                             data['cbm'][0][1] + e_max)
                plt.xlim(0,x_max)
                filename=str(storedir)+str("/")+str('BanddiHSE.png')
                plt.tight_layout()
                plt.savefig(filename,img_format="png")

                plt.close()















                run = Vasprun(ru)
                complete_dos = run.complete_dos
                totup=complete_dos.get_densities(spin=Spin.up)
                totdn=complete_dos.get_densities(spin=Spin.down)
                en = complete_dos.energies
                bandgap=float(complete_dos.get_gap())
                ef = complete_dos.efermi
                en[:] = [x - ef for x in en]



                #allpts = []
                #x=complete_dos.densities
                #y=complete_dos.energies

                plt = get_publication_quality_plot(14, 10)
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.plot(en,totup,'b',linewidth=3)
                plt.plot(en,-totdn,'r',linewidth=3)
                filename=str(storedir)+str("/")+str("TDosHSE.png")
                #print "HSE TDOS file",filename
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename)
                plt.close()

                #print "max den"
                spd_dos = complete_dos.get_spd_dos()
                #print max(spd_dos.densities)
                sdosup=spd_dos[OrbitalType.s].densities[Spin.up]
                pdosup=spd_dos[OrbitalType.p].densities[Spin.up]
                ddosup=spd_dos[OrbitalType.d].densities[Spin.up]
                sdosdn=spd_dos[OrbitalType.s].densities[Spin.down]
                pdosdn=spd_dos[OrbitalType.p].densities[Spin.down]
                ddosdn=spd_dos[OrbitalType.d].densities[Spin.down]
                #try:
                #   fdosup=spd_dos[OrbitalType.f].densities[Spin.up]
                #   fdosdn=spd_dos[OrbitalType.f].densities[Spin.down]
                #except:
                #    pass
                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,sdosup,'r',linewidth=3,label='s')
                plt.plot(en,-sdosdn,'r',linewidth=3)
                plt.plot(en,pdosup,'g',linewidth=3,label='p')
                plt.plot(en,-pdosdn,'g',linewidth=3)
                plt.plot(en,ddosup,'b',linewidth=3,label='d')
                plt.plot(en,-ddosdn,'b',linewidth=3)
                #try:
                #   plt.plot(en,ddosup,'m',linewidth=3,label='f')
                #   plt.plot(en,-ddosdn,'m',linewidth=3)
                #except:
                #   pass
                filename=str(storedir)+str("/")+str("DosHSE.png")
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.legend(prop={'size':26})
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename,pad_inches = 0,bbox_inches='tight')
                plt.close()


                plt = get_publication_quality_plot(14, 10)
                elt_dos = complete_dos.get_element_dos()
                #print "ITEMSSS",elt_dos.items() # = complete_dos.get_element_dos()
                at_symbs=set(run.atomic_symbols)
                #cmap = plt.get_cmap('jet_r')
                for i,sym in enumerate(at_symbs):
                    elt=elt_dos[Element(sym)]
                    #color = cmap(float(i)/len(at_symbs))
                    if i==0:
                        color='b'
                    elif i==1:
                        color='g'
                    elif i==2:
                        color='r'
                    elif i==3:
                        color='c'
                    elif i==4:
                        color='m'
                    elif i==5:
                        color='y'
                    elif i==6:
                        color='k'
                    else:
                        print ("Use different color scheme")
                    plt.plot(en,elt.densities[Spin.up],color=color,linewidth=3,label=sym)
                    plt.plot(en,-elt.densities[Spin.down],color=color,linewidth=3)
                filename=str(storedir)+str("/")+str("EDosHSE.png")
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.legend(prop={'size':26})
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename,bbox_inches='tight')
                plt.close()
                plt = get_publication_quality_plot(14, 10)
                plotter = DosPlotter(zero_at_efermi=True)
                plotter.add_dos_dict(spd_dos)
                #plotter
                filenam=str(storedir)+str("/")+str('Dos1HSE.png')
                #plotter.tight_layout()
                try:
                    plotter.save_plot(filenam,img_format="png",xlim=[-5,10])
                except:
                    plotter.save_plot(filenam,img_format="png")
                #plotter.close()


                spd_dos = complete_dos.get_element_dos()
                plotter = DosPlotter(zero_at_efermi=True)
                plotter.add_dos_dict(spd_dos)
                #plotter
                filenam=str(storedir)+str("/")+str('EDos1HSE.png')
                #plotter.tight_layout()
                ##plt.tight_layout()
                try:
                   plotter.save_plot(filenam,img_format="png",xlim=[-5,10])
                except:
                   plotter.save_plot(filenam,img_format="png")
                #plotter.close()




                filenam=str(storedir)+str("/")+str('BanddiHSE.png')
                bands = run.get_band_structure(kpfile, line_mode = True, efermi = run.efermi)
                if bands.get_band_gap()['direct']==True:
                    bandgap=str(round(bands.get_band_gap()['energy'],3))+str(" D")
                else:
                    bandgap=str(round(bands.get_band_gap()['energy'],3))+str(" I")
                #print ("BANDGAP ISSSS",bands.get_band_gap()['energy'])
#('BANDGAP ISSSS', {'energy': 0.8796999999999997, 'transition': u'K-(0.194,0.194,0.000)', 'direct': False})
                bsp =  BSPlotter(bands)
                #print ("bsp,filename=",bsp,filenam)
                #bsp.save_plot(filenam,smooth=False,img_format="png",ylim=[-5,10])
                plt.close()
                #bsp.close()
                #bsp.save_plot(filename,ylim=(-5,5),img_format="png")
                return bandgap







def banddos(pref='',storedir=None):
    for a in glob.glob(str('*')+str('.json')):
            if a.startswith(pref):
               # from pymatgen.io.vasp.outputs import Vasprun
                from pymatgen.io.vaspio.vasp_output import Vasprun
                ru=str(a.split(".json")[0])+str("/")+str("vasprun.xml")
                kpfile=str(a.split(".json")[0])+str("/")+str("KPOINTS")




    		run = Vasprun(ru, parse_projected_eigen = True)
    		bands = run.get_band_structure(kpfile, line_mode = True, efermi = run.efermi)
    		bsp =  BSPlotter(bands)
    	        zero_to_efermi=True
	        data=bsp.bs_plot_data(zero_to_efermi)
                from pymatgen.util.plotting_utils import get_publication_quality_plot
                plt = get_publication_quality_plot(12, 8)
    #plt = get_publication_quality_plot(12, 8)
                band_linewidth = 3
                x_max = data['distances'][-1][-1]
                #print (x_max)
                for d in range(len(data['distances'])):
                   for i in range(bsp._nb_bands):
                      plt.plot(data['distances'][d],
                             [data['energy'][d]['1'][i][j]
                              for j in range(len(data['distances'][d]))], 'b-',
                             linewidth=band_linewidth)
                      if bsp._bs.is_spin_polarized:
                         plt.plot(data['distances'][d],
                                 [data['energy'][d]['-1'][i][j]
                                  for j in range(len(data['distances'][d]))],
                                 'r--', linewidth=band_linewidth)
                bsp._maketicks(plt)
                if bsp._bs.is_metal():
                     e_min = -10
                     e_max = 10
                     band_linewidth = 3

                for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)

                        for vbm in data['vbm']:
                            plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                        s=100)


        # Main X and Y Labels
                plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
                ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi \
                   else r'$\mathrm{Energy\ (eV)}$'
                plt.ylabel(ylabel, fontsize=30)
                plt.ylim(-5,10)
#    plt.ylim(data['vbm'][0][1] + e_min,
#                             data['cbm'][0][1] + e_max)
                plt.xlim(0,x_max)
                filename=str(storedir)+str("/")+str('Banddi.png')
                plt.tight_layout()
                plt.savefig(filename,img_format="png")

                plt.close()















                run = Vasprun(ru)
                complete_dos = run.complete_dos
                totup=complete_dos.get_densities(spin=Spin.up)
                totdn=complete_dos.get_densities(spin=Spin.down)
                en = complete_dos.energies
                bandgap=float(complete_dos.get_gap())
                ef = complete_dos.efermi
                en[:] = [x - ef for x in en]
                
 
 
                #allpts = [] 
                #x=complete_dos.densities
                #y=complete_dos.energies
             
                plt = get_publication_quality_plot(14, 10)
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.plot(en,totup,'b',linewidth=3)
                plt.plot(en,-totdn,'r',linewidth=3)
                filename=str(storedir)+str("/")+str("TDos.png")
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename)
                plt.close()

                #print "max den"
                spd_dos = complete_dos.get_spd_dos()
                #print max(spd_dos.densities)            
                sdosup=spd_dos[OrbitalType.s].densities[Spin.up] 
                pdosup=spd_dos[OrbitalType.p].densities[Spin.up]
                ddosup=spd_dos[OrbitalType.d].densities[Spin.up]
                sdosdn=spd_dos[OrbitalType.s].densities[Spin.down] 
                pdosdn=spd_dos[OrbitalType.p].densities[Spin.down]
                ddosdn=spd_dos[OrbitalType.d].densities[Spin.down]
                #try: 
                #   fdosup=spd_dos[OrbitalType.f].densities[Spin.up] 
                #   fdosdn=spd_dos[OrbitalType.f].densities[Spin.down] 
                #except:
                #    pass
                plt = get_publication_quality_plot(14, 10)
                plt.plot(en,sdosup,'r',linewidth=3,label='s')
                plt.plot(en,-sdosdn,'r',linewidth=3)
                plt.plot(en,pdosup,'g',linewidth=3,label='p')
                plt.plot(en,-pdosdn,'g',linewidth=3)
                plt.plot(en,ddosup,'b',linewidth=3,label='d')
                plt.plot(en,-ddosdn,'b',linewidth=3)
                #try:
                #   plt.plot(en,ddosup,'m',linewidth=3,label='f')
                #   plt.plot(en,-ddosdn,'m',linewidth=3)
                #except:
                #   pass
                filename=str(storedir)+str("/")+str("Dos.png")
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.legend(prop={'size':26})
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename,pad_inches = 0,bbox_inches='tight')
                plt.close()


                plt = get_publication_quality_plot(14, 10)
                elt_dos = complete_dos.get_element_dos()
                #print "ITEMSSS",elt_dos.items() # = complete_dos.get_element_dos()
                at_symbs=set(run.atomic_symbols)
                #cmap = plt.get_cmap('jet_r')
                for i,sym in enumerate(at_symbs):
                    elt=elt_dos[Element(sym)]
                    #color = cmap(float(i)/len(at_symbs))
                    if i==0:
                        color='b'
                    elif i==1:
                        color='g'
                    elif i==2:
                        color='r'
                    elif i==3:
                        color='c'
                    elif i==4:
                        color='m'
                    elif i==5:
                        color='y'
                    elif i==6:
                        color='k'
                    else:
                        print ("Use different color scheme")
                    plt.plot(en,elt.densities[Spin.up],color=color,linewidth=3,label=sym)
                    plt.plot(en,-elt.densities[Spin.down],color=color,linewidth=3)
                filename=str(storedir)+str("/")+str("EDos.png")
                plt.xlabel("Energies (eV)")
                plt.ylabel("DOS (arb. unit.)")
                plt.legend(prop={'size':26})
                plt.xlim(-5,10)
                plt.tight_layout()
                plt.savefig(filename,bbox_inches='tight')
                plt.close()













                plt = get_publication_quality_plot(14, 10)
                plotter = DosPlotter(zero_at_efermi=True)
                plotter.add_dos_dict(spd_dos)
                #plotter
                filenam=str(storedir)+str("/")+str('Dos1.png')
                #plotter.tight_layout()
                try:
                    plotter.save_plot(filenam,img_format="png",xlim=[-5,10])
                except:
                    plotter.save_plot(filenam,img_format="png")
                #plotter.close()


                spd_dos = complete_dos.get_element_dos()
                plotter = DosPlotter(zero_at_efermi=True)
                plotter.add_dos_dict(spd_dos)
                #plotter
                filenam=str(storedir)+str("/")+str('EDos1.png')
                #plotter.tight_layout()
                ##plt.tight_layout()
                try:
                   plotter.save_plot(filenam,img_format="png",xlim=[-5,10])
                except:
                   plotter.save_plot(filenam,img_format="png")
                #plotter.close()




                filenam=str(storedir)+str("/")+str('Banddi.png')
                bands = run.get_band_structure(kpfile, line_mode = True, efermi = run.efermi)
                if bands.get_band_gap()['direct']==True:
                    bandgap=str(round(bands.get_band_gap()['energy'],3))+str(" D")
                else: 
                    bandgap=str(round(bands.get_band_gap()['energy'],3))+str(" I")
                #print ("BANDGAP ISSSS",bands.get_band_gap()['energy'])
#('BANDGAP ISSSS', {'energy': 0.8796999999999997, 'transition': u'K-(0.194,0.194,0.000)', 'direct': False})
                bsp =  BSPlotter(bands)
                print ("bsp,filename=",bsp,filenam)
                #bsp.save_plot(filenam,smooth=False,img_format="png",ylim=[-5,10])
                plt.close()
                #bsp.close()
                #bsp.save_plot(filename,ylim=(-5,5),img_format="png")
                return bandgap
def elastic_tensor(pref='',storedir=None,name=''):
        ratio_c=1.0
        for a in glob.glob(str('*')+str('.json')):
            if a.startswith("MAIN-RELAX"):
                fold=str(a.split(".json")[0])
                path=str(a.split(".json")[0])+str("/")+str("OUTCAR")
                main_contcar=Structure.from_file(str(a.split(".json")[0])+str("/")+str("POSCAR"))
                old_c=main_contcar.lattice.matrix[2][2]
        for a in glob.glob(str('*')+str('.json')):
            #print "in elastic tensor",a
            if a.startswith(pref):
                fold=str(a.split(".json")[0])
                path=str(a.split(".json")[0])+str("/")+str("OUTCAR")
                contcar=Structure.from_file(str(a.split(".json")[0])+str("/")+str("POSCAR"))
                old_c=contcar.lattice.matrix[2][2]
                filename=str(storedir)+str("/")+str(name)+str('_conv.cif')
                contcar.to(fmt= "cif", filename= filename) 
                

                #ratio_c=0.1*float(abs(contcar.lattice.matrix[2][2]))#*(10**9)*(10**-10) #N/m unit
                try:
                        if "Surf" in a:
                            print "yes",a
                            ratio_c=0.1*float(abs(contcar.lattice.matrix[2][2]))#*(10**9)*(10**-10) #N/m unit
                except:
                        pass
                
                #filename=str(storedir)+str("/")+str("data/")+str(name)+str('_prem.cif')
        print "rato c issssssssssssssssssssssssssssss",ratio_c    
        v=open(path,"r")
        ax_line=''
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
                #print lines[i+4]
                #print lines[i+5]
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
                break
        cwd=str(os.getcwd())
        os.chdir(fold)
        contcar=Structure.from_file("POSCAR")
        symbs=contcar.symbol_set
        tag=str(' '.join([el for el in symbs]))
        mesh=open('meshdos.conf','w')
        if len(contcar.get_primitive_structure())<len(contcar):
           print ("ENFORCINGGG PRIMITIVE AXISSS phoh",SpacegroupAnalyzer(contcar).get_spacegroup_symbol(),len(contcar.get_primitive_structure()),len(contcar))
           finder = SpacegroupAnalyzer(contcar)
           latt=finder.get_lattice_type()
                        ##sgp=str(finder.get_spacegroup_number())
           sgp=str(finder.get_spacegroup_symbol())
           if latt == "rhombohedral":
              transf = np.array([[2, -1, -1], [1, 1, -2], [1, 1, 1]],
                                      dtype=np.float) / 3
           elif "I" in sgp:
                transf = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]],
                              dtype=np.float) / 2

           elif "F" in sgp:
                transf = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                                  dtype=np.float) / 2
           elif "A" in sgp:
              transf = np.array([[2, 0, 0], [0, 1, -1], [0, 1,1 ]],
                                  dtype=np.float) / 2
           elif "C" in sgp:
                if finder.get_crystal_system() == "monoclinic":
                    transf = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
                                      dtype=np.float) / 2
                else:
                    transf = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]],
                                      dtype=np.float) / 2
           else:
                transf = np.eye(3)
           ax=''
           for el in transf:
               ax=str(ax)+str(' ')+str(el[0])+str(' ')+str(el[1])+str(' ')+str(el[2])
           ax_line=str('PRIMITIVE_AXIS = ')+str(ax)+'\n'
           #try:
           #   mesh.write(ax_line)
           #except:
           #   pass
        line=str('FORCE_CONSTANTS = READ')+'\n'
        mesh.write(line)
        line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
        mesh.write(line)
        line=str('ATOM_NAME = ')+str(tag)+'\n'
        mesh.write(line)
        line=str('DIM = 1 1 1')+'\n'
        mesh.write(line)
        line=str('MP = 31 31 31')+'\n'
        mesh.write(line)
        mesh.close()





        incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = 45678,
           NBANDS=456789,
           ISMEAR = 0,
           EDIFF = '1E-7',
           LCHARG = '.FALSE.',
           NEDOS = 5000,
           ISPIN = 2,
           ISIF = 2,
           IBRION = 1,
           NELM = 400,
           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
        incar = Incar.from_dict(incar_dict)
        user_incar_settings={"EDIFF":1E-6,"ISIF":2,"NSW":0,"LORBIT":11,"ENCUT":34,"LWAVE":'.FALSE.',"PREC":'Accurate'}
        mpvis = MPNonSCFVaspInputSet(user_incar_settings=incar)
        kpoints=mpvis.get_kpoints(contcar)
        all_kp=kpoints.kpts
        labels=kpoints.labels
        all_labels=''
        allline=''
        for l in labels:
           if l=='':
              l=None

           all_labels=all_labels+str(l)+str(' ')
        for k in all_kp:
           #allline=allline+str(k)
           allline=allline+ str(k[0])+str(' ')+str(k[1])+str(' ')+str(k[2])+str(' ')
        #print allline
        file=open('bandd.conf','w')
        line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
        file.write(line)
        if len(contcar.get_primitive_structure())<len(contcar):
           try:
              file.write(ax_line)
           except:
               pass
        line=str('ATOM_NAME = ')+str(tag)+'\n'
        file.write(line)
        line=str('DIM = 1 1 1')+'\n'
        file.write(line)
        line=str('FORCE_CONSTANTS = READ')+'\n'
        file.write(line)
        line=str("BAND= ")+str(allline)+'\n'
        file.write(line)
        line=str("BAND_LABELS= ")+str(all_labels)+'\n'
        file.write(line)
        file.close()




        pg_phonopy='na'


#check if there is mesh.yaml
        meshf=str(os.getcwd())+str("/")+str("mesh.yaml")
        forcf=str(os.getcwd())+str("/")+str("FORCE_CONSTANTS")
        irf=str(os.getcwd())+str("/")+str("irreps.yaml")
        thermf=str(os.getcwd())+str("/")+str("thermal_properties.yaml")
        bandf=str(os.getcwd())+str("/")+str("band.yaml")

        #if not os.path.exists(meshf):
        try:
           if not os.path.exists(forcf):
              cmd=str('phonopy --fc vasprun.xml')
              os.system(cmd)
        except:
           pass


        try:
           cmd=str('rm irreps.yaml')
           #os.system(cmd)
        except:
            pass
        try:

           if not os.path.exists(meshf):
              cmd=str('phonopy -p meshdos.conf')
              os.system(cmd)
              cmd=str('phonopy -t meshdos.conf')
              os.system(cmd)
           if not os.path.exists(bandf):
              cmd=str('phonopy -p bandd.conf')
              os.system(cmd)
          # if not os.path.exists(irf):
          #    cmd=str('phonopy -p mesh.conf')
          #    os.system(cmd)
        except:
            pass



        try:
           cmd=str('phonopy -p meshdos.conf')
           #os.system(cmd)
           cmd=str('phonopy -t meshdos.conf')
           #os.system(cmd)
           cmd=str('phonopy -p bandd.conf')
           #os.system(cmd)
           cmd=str('phonopy --symmetry >PHONOPY_SYMMETRY')
           os.system(cmd)
           fp=open("PHONOPY_SYMMETRY","r")
           for lp in fp:
               if "point_group_type" in lp:
                   pg_phonopy=str(lp) #.split('point_group_type:')[1]
                   print ("pg_phonopy",pg_phonopy)
           fp.close()#=open("PHONOPY_SYMMETRY","r")
        except:    
           pass


        filename=str(storedir)+str("/")+str("PBAND.png")
        try:
          cmd=str('bandplot -o ')+str(filename)
          os.system(cmd)
          #cmd=str('phonopy --symmetry >PHONOPY_SYMMETRY')
          #os.system(cmd)
        except:
           pass


        #fp=open("PHONOPY_SYMMETRY","r")
        #for lp in fp:
        #    if "point_group_type" in lp:
        #        print line
        #        pg_phonopy=str(lp)#.split('point_group_type:')[1]
        #        #pg_phonopy=line.split('point_group_type:')[1]
        #        print ("pg_phonopy",pg_phonopy)
        #fp.close()#=open("PHONOPY_SYMMETRY","r")


        try:
		with open('mesh.yaml', 'r') as f:
		   doc = yaml.load(f)
		nmodes=doc['phonon'][0]['band']
		ph_modes=[]
		
		for p in nmodes:
		    ph_modes.append(p['frequency'])
		ph_modes=sorted(set(ph_modes))
		f=open('total_dos.dat','r')
		freq=[]
		pdos=[]
		for lines in f.readlines():
		    if str(lines.split()[0]).startswith("#"):
		       print lines
		    else:
		       freq.append(float(lines.split()[0]))
		       #freq.append(float(33.35641)*float(lines.split()[0]))
		       pdos.append(float(lines.split()[1]))
	    #plt.xlim(0,300)
		#plt.close()
		plt = get_publication_quality_plot(14, 10)
		plt.xlabel("Frequency ($cm^{-1}$)")
		plt.ylabel("PDOS (arb. unit.)")
		plt.plot(freq,pdos,'b',linewidth=2)
		filename=str(storedir)+str("/")+str("pdos_conv.png")
		plt.tight_layout()
		plt.savefig(filename)
		plt.close()
        except:
           pass
        os.chdir(cwd)



        return c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,ph_modes,pg_phonopy

def  write_html(storedir='',name='',nmodes=None,c11=0,c12=0,c13=0,c14=0,c15=0,c16=0,c21=0,c22=0,c23=0,c24=0,c25=0,c26=0,c31=0,c32=0,c33=0,c34=0,c35=0,c36=0,c41=0,c42=0,c43=0,c44=0,c45=0,c46=0,c51=0,c52=0,c53=0,c54=0,c55=0,c56=0,c61=0,c62=0,c63=0,c64=0,c65=0,c66=0):

	f=open(str(storedir)+str("/")+str(name)+str(".html"),"w")
	line=str("<!DOCTYPE html>")+'\n'
	f.write(line)
	line=str("<html>")+'\n'
	f.write(line)
	line=str("<title>JSmol-MathJax Compatibility Test</title>")+'\n'
	f.write(line)
	line=str("<head>")+'\n'
	f.write(line)
	line=str("</script>")+'\n'
	f.write(line)
	line=str('<script type="text/javascript"  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>')+'\n'
	f.write(line)
	line=str('<script type="text/javascript" src="JSmol.min.js"></script>')+'\n'
	f.write(line)
	line=str('<script type="text/javascript">')+'\n'
	f.write(line)
	line=str("Info = {")+'\n'
	f.write(line)
	line=str("	width: 400,")+'\n'
	f.write(line)
	line=str("	height: 400,")+'\n'
	f.write(line)
	line=str("	debug: false,")+'\n'
	f.write(line)
	line=str('	color: "0xC0C0C0",')+'\n'
	f.write(line)
	line=str("  disableJ2SLoadMonitor: false,")+'\n'
	f.write(line)
	line=str("  disableInitialConsole: true,")+'\n'
	f.write(line)
	line=str("	addSelectionOptions: false,")+'\n'
	f.write(line)
	line=str('	use: "HTML5",')+'\n'
	f.write(line)
	line=str("	readyFunction: null,")+'\n'
	f.write(line)
	#line=str('	script: "load data/')+str(name)+str('_prem.cif"')+'\n'
	line=str('	script: "load VASP-FIGS/')+str(name)+str("/")+str(name)+str('_prem.cif"')+'\n'
	f.write(line)
	line=str("}")+'\n'
	f.write(line)
	line=str("loadJmol = function() {")+'\n'
	f.write(line)
	line=str('  $("#mathdiv").show();')+'\n'
	f.write(line)
	line=str('	$("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))')+'\n'
	f.write(line)
	line=str("}")+'\n'
	f.write(line)
	line=str("$(document).ready(function(){")+'\n'
	f.write(line)
	line=str("  checkMathJax(loadJmol);")+'\n'
	f.write(line)
	line=str("});")+'\n'
	f.write(line)
	line=str("checkMathJax = function(fReady) {")+'\n'
	f.write(line)
	line=str('  return (!$("#MathJax_Message")[0] ')+'\n'
	f.write(line)
	line=str("    || $(")+str("'")+str('#MathJax_Message:contains("Loading")')+str("'")+str(")[0]")+'\n'
	f.write(line)
	line=str("    ? setTimeout(function(){checkMathJax(fReady)},10)")+'\n'
	f.write(line)
	line=str("    : fReady());")+'\n'
	f.write(line)
	line=str("}")+'\n'
	f.write(line)
	line=str("</script>")+'\n'
	f.write(line)
	line=str("</head>")+'\n'
	f.write(line)
	line=str("<body>")+'\n'
	f.write(line)
	line=str('<table><tr><td style="width:350px">')+'\n'
	f.write(line)
	#line=str("</td><td>")+'\n'
	#f.write(line)
	line=str("<span id=mydiv></span>")+'\n'
	f.write(line)
	line=str("</td></tr></table>")+'\n'
	f.write(line)






	line=str('<h2 style="color:blue;">Convergence</h2>')+'\n'
	f.write(line)
        line=str("<p> This is convergenece test </p>")+'\n'
	f.write(line)
	line=str('    <div class="row">')+'\n'
	f.write(line)
	line=str('      <div class="col s12 m6">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/Encut.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('      </div>')+'\n'
	f.write(line)


	line=str("<h2>Structural analysis</h2>")+'\n'
	f.write(line)
	line=str('    <div class="row">')+'\n'
	f.write(line)
	line=str('      <div class="col s12 m6">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/xrd.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/rdf.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('      </div>')+'\n'
	f.write(line)




	line=str('<h2>Electronic structure</h2>')+'\n'
	f.write(line)
	line=str('      <div class="col s12 m6">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/TDos.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/Banddi.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('      </div>')+'\n'
	f.write(line)
	line=str("<h2>Optical properties</h2>")+'\n'
	f.write(line)
	line=str('      <div class="col s12 m6">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/Real.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/Imag.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str("      </div>")+'\n'
	f.write(line)
	line=str("<h2>Phonon properties</h2>")+'\n'
	f.write(line)
	line=str('      <div class="col s12 m6">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/pdos_prim.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('        <image src="VASP-FIGS/')+str(name)+str('/pdos_conv.png" width="450" height="350">')+'\n'
	f.write(line)
	line=str('      </div>')+'\n'
	f.write(line)
	line=str("    </div>")+'\n'
	f.write(line)
	line=str("<h2>Phonon</h2>")+'\n'
	f.write(line)
	line=str("<table>")+'\n'
	f.write(line)
	line=str("  <tr>")+'\n'
	f.write(line)
	line=str("    <th>Mode</th>")+'\n'
	f.write(line)
	line=str("    <th>Name</th>")+'\n'
	f.write(line)
        for modes in nmodes:
            try:
        	line=str("                <tr>")+'\n'
	        f.write(line)
	        line=str("                    <td>")+str(modes['frequency'])+str("</td>")+'\n'
        	f.write(line)
	        line=str("                    <td>")+str(modes['ir_label'])+str("</td>")+'\n'
        	f.write(line)
        	line=str("                </tr>")+'\n'
	        f.write(line)
            except:
        	line=str("                <tr>")+'\n'
	        f.write(line)
	        line=str("                    <td>")+str(modes['frequency'])+str("</td>")+'\n'
        	f.write(line)
        	line=str("                </tr>")+'\n'
	        f.write(line)
	line=str("  </tr>")+'\n'
	f.write(line)
	line=str("</table>")+'\n'
	f.write(line)


	line=str("</body>")+'\n'
	f.write(line)
	line=str("</body>")+'\n'
	f.write(line)





	f.close()
def Ramanplot(file='',storedir=''):
    plt = get_publication_quality_plot(14, 10)
    f=open(file,'r')
    freq=[]
    raman=[]
    for lines in f.readlines():
        if str(lines.split()[0]).startswith("#"):
           print lines
        else:
           freq.append(float(lines.split()[1]))
           raman.append(float(lines.split()[4]))
    #plt.xlim(50,300)
    plt.xlabel("Wavenumber($cm^{-1}$)")
    plt.ylabel("Raman Intensity")
#plt.plot(freq,raman,'b',linewidth=2)
    plt.bar(freq,raman,width=5,color='b')
    filename=str(storedir)+str("/")+str("Raman.png")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    f.close()


def phonon(storedir=None,name=None):
    cwd=str(os.getcwd())
    print ("phonn cwd=",cwd)
    for a in glob.glob(str('*')+str('.json')):
            #print a
            if a.startswith('MAIN-ELAST'):
                fold=str(a.split(".json")[0])
                path=str(a.split(".json")[0])+str("/")+str("OUTCAR")
        #ramandat=str(os.getcwd())+str('/')+a+str('/vasp_raman.dat')
        #print ("RAMiANDAT====",ramandat,storedir)
    try:
       Ramanplot(file=ramandat,storedir=storedir)
    except:
       pass
    phdir=str(os.getcwd())+str('/')+str(fold)
    #print ("phdir=",phdir)
    os.chdir(phdir)
    #print ("RAMANDAT2====",ramandat,phdir)
    contcar=Structure.from_file("POSCAR")


    symbs=contcar.symbol_set
    tag=str(' '.join([el for el in symbs]))
    mesh=open('meshdos.conf','w')
    line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
    mesh.write(line)
    line=str('FORCE_CONSTANTS = READ')+'\n'
    mesh.write(line)
    try:
       cmd=str('rm irreps.yaml')
       #os.system(cmd)
    except:
        pass
    ax_line=str('')
    if len(contcar.get_primitive_structure())<len(contcar):
       print ("ENFORCINGGG PRIMITIVE AXISSSpho222",SpacegroupAnalyzer(contcar).get_spacegroup_symbol())
       finder = SpacegroupAnalyzer(contcar)
       latt=finder.get_lattice_type()
                        ##sgp=str(finder.get_spacegroup_number())
       sgp=str(finder.get_spacegroup_symbol())
       if latt == "rhombohedral":
              transf = np.array([[2, -1, -1], [1, 1, -2], [1, 1, 1]],
                                      dtype=np.float) / 3
       elif "I" in sgp:
            transf = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]],
                              dtype=np.float) / 2

       elif "F" in sgp:
            transf = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                              dtype=np.float) / 2
       elif "A" in sgp:
              transf = np.array([[2, 0, 0], [0, 1, -1], [0, 1,1 ]],
                                  dtype=np.float) / 2
       elif "C" in sgp:
            if finder.get_crystal_system() == "monoclinic":
                transf = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
                                  dtype=np.float) / 2
            else:
                transf = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]],
                                  dtype=np.float) / 2
       else:
            transf = np.eye(3)
       ax=''
       for el in transf:
           ax=str(ax)+str(' ')+str(el[0])+str(' ')+str(el[1])+str(' ')+str(el[2])
       ax_line=str('PRIMITIVE_AXIS = ')+str(ax)+'\n'
       try:
          mesh.write(ax_line)
       except:
          pass
    line=str('ATOM_NAME = ')+str(tag)+'\n'
    mesh.write(line)
    line=str('DIM = 1 1 1')+'\n'
    mesh.write(line)
    line=str('MP = 31 31 31')+'\n'
    mesh.write(line)
    mesh.close()



    anim=open('anim.conf','w')
    line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
    anim.write(line)
    line=str('FORCE_CONSTANTS = READ')+'\n'
    anim.write(line)
    line=str('ATOM_NAME = ')+str(tag)+'\n'
    anim.write(line)
    line=str('ANIME_TYPE = JMOL')+'\n'
    anim.write(line)
    line=str('ANIME = 0 7 0  0 0 0')+'\n'
    anim.write(line)
    line=str('DIM = 1 1 1')+'\n'
    anim.write(line)
    #line=str('MP = 31 31 31')+'\n'
    #anim.write(line)
    anim.close()


    

    vfile=str(storedir)+str("/")+str(name)+str("_vphonon.html")
    vp=open(vfile,"w")
    line=str('<!DOCTYPE html>')+'\n'
    vp.write(line)
    line=str('<html>')+'\n'
    vp.write(line)
    line=str('  <head>')+'\n'
    vp.write(line)
    line=str('    <meta charset="utf-8">')+'\n'
    vp.write(line)
    line=str('    <title>VPHONON</title>')+'\n'
    vp.write(line)
    line=str('    <p> Right click on the atomic structure. Then select third option "model" and then select the frequency to be visualized.</p>')+'\n'
    vp.write(line)
    line=str('    <script src="JSmol.min.js"></script>')+'\n'
    vp.write(line)
    line=str('    <script src="js/Jmol2.js"></script>')+'\n'
    vp.write(line)
    line=str('    <script>')+'\n'
    vp.write(line)
    line=str('//jmolInitialize(".");')+'\n'
    vp.write(line)
    line=str('function loadVibs(){')+'\n'
    vp.write(line)
    line=str(' var Info = jmolGetPropertyAsArray("auxiliaryinfo.VibFreqs")')+'\n'
    vp.write(line)
    line=str(' if(!Info){')+'\n'
    vp.write(line)
    line=str('	alert("No VibFreqs")')+'\n'
    vp.write(line)
    line=str('	return')+'\n'
    vp.write(line)
    line=str(' }')+'\n'
    vp.write(line)
    line=str(' var s="<select id=vib onchange=')+str("'showFrame(value)' onkeypress='showFrame()'><option value='1'>")+str('no vibration (model 1)</option>";')+'\n'
    vp.write(line)
    line=str(' for(var i=0;i<Info.length;i++)')+'\n'
    vp.write(line)
    line=str('  s+="<option value='"+(i+2)+"'>"+Info[i].freq + " cm-1 "+Info[i].label+"</option>"')+'\n'
    vp.write(line)
    line=str(' s+="</select>"')+'\n'
    vp.write(line)
    line=str(' document.getElementById("freqDiv").innerHTML = s')+'\n'
    vp.write(line)
    line=str('}')+'\n'
    vp.write(line)
    line=str('function showFrame(i){')+'\n'
    vp.write(line)
    line=str(' if(arguments.length == 0){')+'\n'
    vp.write(line)
    line=str('    setTimeout("showFrame(-1)",100)')+'\n'
    vp.write(line)

    line=str('    return')+'\n'
    vp.write(line)
    line=str(' }')+'\n'
    vp.write(line)
    line=str(' if(i==-1) {')+'\n'
    vp.write(line)
    line=str('  var d=document.getElementById("vib")')+'\n'
    vp.write(line)
    line=str('  i=d[d.selectedIndex].value')+'\n'
    vp.write(line)
    line=str(' }')+'\n'
    vp.write(line)
    line=str(' jmolScript("frame " + i)')+'\n'
    vp.write(line)
    line=str('}')+'\n'
    vp.write(line)
    line=str('    </script>')+'\n'
    vp.write(line)
    line=str('  </head>')+'\n'
    vp.write(line)
    line=str('  <body>')+'\n'
    vp.write(line)
    line=str('    <form>')+'\n'
    vp.write(line)
    line=str('      <table border="1">')+'\n'
    vp.write(line)
    line=str('	<tr>')+'\n'
    vp.write(line)
    line=str('	  <td>')+'\n'
    vp.write(line)
    line=str('	    <script>')+'\n'
    vp.write(line)
    line=str('	      var script =')+'\n'
    vp.write(line)
    line=str('	      "load ')+str('VASP-FIGS/')+str(name)+str('/anime.xyz_jmol; wireframe .1; cpk off;" +')+'\n'
    #line=str('	      "load data/anime.xyz_jmol; wireframe .1; cpk off;" +')+'\n'
    vp.write(line)
    line=str('	      "frame 1; vectors 0.1;vector scale 5.0 ; color vectors yellow;" +')+'\n'
    vp.write(line)
    line=str(' 	      "move 10 -20 10 0 0 0 0 0 1; delay 1; vibration on;vibration scale 0.4;";')+'\n'
    vp.write(line)
    line=str(' 	      jmolApplet(500, script);')+'\n'
    vp.write(line)
    line=str('  </script>')+'\n'
    vp.write(line)
    line=str('	  </td>')+'\n'
    vp.write(line)
    line=str(' 	  <td>')+'\n'
    vp.write(line)
    line=str(' 	    </blockquote>')+'\n'
    vp.write(line)
    line=str(' 	    <blockquote>')+'\n'
    vp.write(line)
    line=str(' 	      <script>')+'\n'
    vp.write(line)
    line=str(' 		jmolCheckbox("vibration on", "vibration off", "vibration", "checked");')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 		jmolCheckbox("vectors 0.05", "vectors off", "vectors", "checked");')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 		jmolRadioGroup([')+'\n'
    vp.write(line)
    line=str(' 		["color vectors yellow", null, "checked"],')+'\n'
    vp.write(line)
    line=str(' 		"color vectors purple"')+'\n'
    vp.write(line)
    line=str(' 		]);')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 		jmolRadioGroup([')+'\n'
    vp.write(line)
    line=str(' 		["spacefill off", null, "checked"],')+'\n'
    vp.write(line)
    line=str(' 		"spacefill 20%",')+'\n'
    vp.write(line)
    line=str(' 		"spacefill 100%"')+'\n'
    vp.write(line)
    line=str(' 		]);')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 		jmolRadioGroup([')+'\n'
    vp.write(line)
    line=str(' 		"wireframe on",')+'\n'
    vp.write(line)
    line=str(' 		["wireframe 0.1", null, "checked"],')+'\n'
    vp.write(line)
    line=str(' 		]);')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 		jmolCheckbox("spin on", "spin off", "spin");')+'\n'
    vp.write(line)
    line=str(' 		jmolBr();')+'\n'
    vp.write(line)
    line=str(' 	      </script>')+'\n'
    vp.write(line)
    line=str(' 	    </blockquote>')+'\n'
    vp.write(line)
    line=str(' 	    </td>')+'\n'
    vp.write(line)
    line=str(' 	  </tr>')+'\n'
    vp.write(line)
    line=str('      </table>')+'\n'
    vp.write(line)
    line=str('    </form>')+'\n'
    vp.write(line)
    line=str('</html>')+'\n'
    vp.write(line)
    vp.close()
    #line=str(' }')+'\n'
    #vp.write(line)










    mesh=open('mesh.conf','w')
    line=str('IRREPS = 0 0 0 1e-3')+'\n'
    mesh.write(line)
    mesh.write(ax_line)
    try:
        mesh.write(ax_line)
    except:
         pass
    line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
    mesh.write(line)
    line=str('SHOW_IRREPS = .TRUE.')+'\n'
    mesh.write(line)
    line=str('FORCE_CONSTANTS = READ')+'\n'
    mesh.write(line)
    line=str('DIM = 1 1 1')+'\n'
    mesh.write(line)
    line=str('MP = 31 31 31')+'\n'
    mesh.write(line)
    mesh.close()




#check if there is mesh.yaml
    forcf=str(os.getcwd())+str("/")+str("FORCE_CONSTANTS")
    meshf=str(os.getcwd())+str("/")+str("mesh.yaml")
    irf=str(os.getcwd())+str("/")+str("irreps.yaml")
    thermf=str(os.getcwd())+str("/")+str("thermal_properties.yaml")
    bandf=str(os.getcwd())+str("/")+str("band.yaml")

    #if not os.path.exists(meshf):
    try:
       if not os.path.exists(forcf):
         cmd=str('phonopy --fc vasprun.xml')
         os.system(cmd)
    except:
       pass


    try:
       cmd=str('rm irreps.yaml')
       #os.system(cmd)
    except:
        pass
    try:

       if not os.path.exists(meshf):
          cmd=str('phonopy -p meshdos.conf')
          os.system(cmd)
          cmd=str('phonopy -t meshdos.conf')
          os.system(cmd)
       if not os.path.exists(bandf):
          cmd=str('phonopy -p bandd.conf')
          os.system(cmd)
       if not os.path.exists(irf):
          cmd=str('phonopy -p mesh.conf --tolerance=1e-3')
          os.system(cmd)
    except:
        pass




    try:
       cmd=str('phonopy -p mesh.conf --tolerance=1e-3')
       #os.system(cmd)
    except:
       pass
    try:
       cmd=str('phonopy -p meshdos.conf')
       #os.system(cmd)
    except:
        pass
    try:
       cmd=str('phonopy -p anim.conf')
       os.system(cmd)
       filename1=str(os.getcwd())+str("/")+str("anime.xyz_jmol")
       filename2=str(storedir)+str("/")+str(name)+str("/")+str("anime.xyz_jmol")
       #print (filename1)
       #print (filename2)
       shutil.copy2(filename1,filename2)
    except:
        pass
    try:
       with open('irreps.yaml', 'r') as f:
           doc = yaml.load(f)
       nmodes=doc['normal_modes']
    except:
       pass
          
    for modes in nmodes:
        try:
           print modes['band_indices'],modes['frequency'],modes['ir_label']
        except:
           pass
    f=open('total_dos.dat','r')
    freq=[]
    pdos=[]
    for lines in f.readlines():
        if str(lines.split()[0]).startswith("#"):
           print lines
        else:
           freq.append(float(33.35641)*float(lines.split()[0]))
           pdos.append(float(lines.split()[1]))
    #plt.xlim(0,300)
    #plt.close()
    plt = get_publication_quality_plot(14, 10)
    plt.xlabel("Frequency ($cm^{-1}$)")
    plt.ylabel("Phonon DOS (arb. unit.)")
    plt.plot(freq,pdos,'b',linewidth=2)
    plt.tight_layout()
    filename=str(storedir)+str("/")+str("pdos_prim.png")
    plt.savefig(filename)
    #plt.close()

    plt.close()
    os.chdir(cwd)
    print os.getcwd()




    return nmodes


def main(chdir=None,storedir='/data/knc6/NIST2/OLD/STORE'):
    os.chdir(chdir)
    try:
        K_fig(pref='DENbulk',storedir=storedir)
    except:
      pass
    try:
        EN_fig(pref='ENCUTbulk',storedir=storedir)
    except:
        pass
    try:

        banddos(pref='MAIN-BAND',storedir=storedir)
    except:
        pass
    try:
        c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66=elastic_tensor(pref='MAIN-ELASTIC',storedir=storedir)
    except:
        pass
    try:
       nmodes=phonon(storedir=storedir)
    except:
      pass
    try:
       boltz_bader(pref='MAIN-RELAX',storedir=storedir)
    except:
       pass
    os.chdir('../')

def loop_over():
   count=0
   key=0
   df=open("DFT_dathold32g.csv","w")
   line=str('key,')+str('symbol,')+str('period,')+str('group')+'\n'
   Infarray=[]
   df.write(line)
   new_cal=[]
   all_folds=[]
   for file in glob.glob("*-*/*_*"):

       if os.path.isdir(file):
           all_folds.append(file)
   hold=sorted( all_folds, key=os.path.getmtime)
   #hold=['Solar-SemiRaritan20/mvc-2597_PBEBO'] #sorted( all_folds, key=os.path.getmtime)
   #hold=['TE-bulk/mp-1201_bulk_PBEBO'] #sorted( all_folds, key=os.path.getmtime)
   #hold=['2D-1L/POSCAR-mp-558544-1L.vasp_PBEBO'] #sorted( all_folds, key=os.path.getmtime)
   #hold=['2D-1L/POSCAR-mp-984-1L_PBEBO'] #sorted( all_folds, key=os.path.getmtime)
   #hold=['2D-1L/POSCAR-mp-1634-1L.vasp_PBEBO'] #sorted( all_folds, key=os.path.getmtime)
   #hold=['Elements-bulkk/mp-149_bulk_PBE','2D-1L/POSCAR-mp-1634-1L.vasp_PBEBO','2D-bulk/mp-66_bulk_LDA'] #sorted( all_folds, key=os.path.getmtime)
   print (hold),len(hold)
   xc_info_arr=[]
   #for file in glob.glob( 'MoS2_mp-1434_PBE'):
   #for file in glob.glob("/data/knc6/NIST4/Elements-bulkk/mp-149_bulk_LDA"):
   #for file in glob.glob("2D-1L/*_*"):
   #for file in glob.glob("Elements-bulkk/mp-*"):
   #for file in glob.glob("2D-*/*_*"):
   #for file in glob.glob("*-*/*_*"):
   #for file in glob.glob("/data/knc6/NIST4/2D-bulk/mp-66_bulk_LDA"):
   #for file in glob.glob("/data/knc6/NIST4/Elements-bulkk/mp-149_bulk_PBE"):
   for tmp_file in hold:
        all_refs=[]
        xc_info={}
        cal=file.split('_')[0]
        file=str(os.getcwd())+str("/")+str(tmp_file)
        cwd=str(os.getcwd())
        os.chdir(file)
        try:
          cmmd=str('rm *json.json*')
          os.system(cmmd)
        except:
          pass
        print ("FILE CWDDDD",file,os.getcwd())   
        done=False
        partial_done=False
        info=[]
        dirs_list=[d for d in os.listdir('./') if os.path.isdir(os.path.join('./', d))]
        for l in dirs_list:
        #for ub_file in glob.glob('MAIN-*'):
        #    if "MAIN-BAND" in ub_file:
            if "MAIN-BAND" in l:
               fb_name=str(l)+str(".json")
               if not os.path.isfile(fb_name):
                  fb=open(fb_name,'w')
                  fb.close()
            if "MAIN-RELAX" in l:
                partial_done=True
                print "partial done"
                p_strt=Structure.from_file(str(l.split(".json")[0])+str("/")+str("POSCAR"))
                fb_name=str(l)+str(".json")
                if not os.path.isfile(fb_name):
                   fb=open(fb_name,'w')
                   fb.close()
        for sub_file in glob.glob('*.json'):
            #print sub_file
        #for sub_file in glob.glob(str(file)+str('/*.json')):
            if 'MAIN-BAND' in sub_file:
               #print ("MAINBAND FILE=",(str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml")))
               #done=(Vasprun(str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"))).converged
               #print "MAIN-BAND isssss",sub_file
               #print (Vasprun(str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"), parse_projected_eigen=False)).converged
            #if 'MAIN-ELASTIC' in sub_file:
            #if sub_file.startswith('MAIN-RELAX'):
               try:
               #if  os.path.isfile((str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"))):
                   if  os.path.isfile((str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"))):
                       done=(Vasprun(str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"))).converged
                       #outcar=str(os.getcwd())+str("/")+str(sub_file.split(".json")[0])+str("/")+str("OUTCAR")
                       #fout=open(outcar,'r')
                       #print "outcar isss",outcar
                       #for fl in fout:
                       #   if   "General timing and accounting informations for this job" in fl :
                       #        done=True
                       #        print ("Proceedding ahead")
                       #fout.close()
                       strt=Structure.from_file(str(sub_file.split(".json")[0])+str("/")+str("POSCAR"))
                       red_formula=strt.composition.reduced_formula
                       try:
                         final_enp='na'
                         final_enp=round(float((Vasprun(str(sub_file.split(".json")[0])+str("/")+str("vasprun.xml"))).final_energy)/float(strt.composition.num_atoms),4)
                         print ("FINAL ENNNNP ISSSSSSSSSSSSSSSSSSSSSSLLSLSLSLSLSLS=",final_enp)
                       except:
                          pass
                           
                       try:
                           finder = SpacegroupAnalyzer(strt)
                        ##sgp=str(finder.get_spacegroup_number())
                           sgp=str(finder.get_spacegroup_symbol())
                           abc=str('na')
                           angles=str('na')
                           abc=str("a ")+ str(round(strt.lattice.abc[0],3))+str(" &#8491; , b ")+str(round(strt.lattice.abc[1],3))+str(" &#8491; , c ")+str(round(strt.lattice.abc[2],3))+str(" &#8491;")
                           angles=str("&alpha; ")+ str(round(strt.lattice.angles[0],3))+str(" &#176; , &beta; ")+str(round(strt.lattice.angles[1],3))+str(" &#176; , &gamma; ")+str(round(strt.lattice.angles[2],3))+str(" &#176;")
                       except:
                          pass
               except:
                 pass
        if done==True :
                        print ("done true",os.getcwd())
                        formula=strt.composition.formula
                        red_formula=strt.composition.reduced_formula
                        print "formula=",formula
			new_cal.append(file.split('_')[0])
			for cc in new_cal:
			   ccc=0
			   if cal in new_cal:
			      ccc=ccc+1

		  # for file in glob.glob('*_mp*'):
			group=0
			period=ccc

			#print (file)
			count=count+1
			key=key+1
			#storedir=str('/data/knc6/NIST2/OLD/STORE')
                        with open('JARVIS-ID','r') as jf:
                             jl=jf.read()
                             name=jl.split('\n')[0] #str("JVASP-")+str(count)
			storedir=str('/rk2/knc6/FRESH32g_DIR')
			#storedir=str('/data/knc6/STORAGE/NEW_STORAGE')
			dir_name=os.path.join(storedir,str(name))
                        if not os.path.exists(dir_name): os.makedirs(dir_name)
			##take care of cwd fin_en,search,num,avg,additional_info=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
                        #print "Effective mass",avg['p'][300][0]
                        #print "Effective mass",avg['p'][300][0][0][0],avg['p'][300][0][0][1],avg['p'][300][0][0][2]
                        #print "Effective mass",avg['p'][300][0][1][0],avg['p'][300][0][1][1],avg['p'][300][0][1][2]
                        #print "Effective mass",avg['p'][300][0][2][0],avg['p'][300][0][2][1],avg['p'][300][0][2][2]
			print "COUNTTTT===",count,dir_name
                        ref_numb='na'
			#vtot(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
			#fin_en,search=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JV2D-")+str(count))
			#if  os.path.exists(cwd+'/'+str(file)+str("/")+str("REFERENCE")): 
                        #   ref=open(cwd+'/'+str(file)+str("/")+str("TAGS")+str("REFERENCE"),"r")
                        #   for line in ref:
                        #       ref_numb=(str(line)).split('\n')[0]
                        #   ref.close()
			   #shutil.copy2(cwd+'/'+str(file)+str("/")+str("REFERENCE"),dir_name)
			tag='na'
                        print "tag",cwd+str('/')+str(file)+str("/")+str("FUNCTIONAL")
			if  os.path.exists("FUNCTIONAL"):
			#if  os.path.exists(cwd+str('/')+str(file)+str("/")+str("FUNCTIONAL")):
			    tag_file= open("FUNCTIONAL","r")
			    for line in tag_file:
				tag=(str(line)).split('\n')[0]
                                print ("tag====",tag)
			    tag_file.close()
			if tag !='na':
			   print ("New tag assigned")
			else:
			  tag=str("Unknown")
                        print "TAGGGG IS ",tag
			#name=str("JVASP-")+str(count)
			f=open(str(storedir)+str("/")+str(name)+str(".html"),"w")
			line=str("<!DOCTYPE html>")+'\n'
			f.write(line)
			line=str("<html>")+'\n'
			f.write(line)
			line=str('<body style="background-color:lightgreen;">')+'\n'
			f.write(line)
			line=str("<title>JARVIS-DFT</title>")+'\n'
			#line=str("<title>JSmol-MathJax Compatibility Test</title>")+'\n'
			f.write(line)



			line=str('<style>')+'\n'
			f.write(line)
			line=str('ul {')+'\n'
			f.write(line)
			line=str('    list-style-type: none;')+'\n'
			f.write(line)
			line=str('    margin: 0;')+'\n'
			f.write(line)
			line=str('    padding: 0;')+'\n'
			f.write(line)
			line=str('    overflow: hidden;')+'\n'
			f.write(line)
			line=str('    background-color: 	#808080;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)

			line=str('li {')+'\n'
			f.write(line)
			line=str('    float: left;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)

			line=str('li a {')+'\n'
			f.write(line)
			line=str('    display: inline-block;')+'\n'
			f.write(line)
			line=str('    color: white;')+'\n'
			f.write(line)
			line=str('    text-align: center;')+'\n'
			f.write(line)
			line=str('    padding: 14px 16px;')+'\n'
			f.write(line)
			line=str('    text-decoration: none;')+'\n'
			f.write(line)
			line=str('}')+'\n'

			line=str('li a:hover {')+'\n'
			f.write(line)
			line=str('    background-color:#000;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)
			line=str('</style>')+'\n'
			f.write(line)

			line=str('<ul>')+'\n'
			f.write(line)

				
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/JARVIS.html">HOME</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/JVASP.html">JARVIS-DFT</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/periodic.html">JARVIF-FF</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/documentation.html">DOCUMENTATION</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/others.html">OTHER RESOURCES</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/contact.html">CONTACT</a></li>')+'\n'
			f.write(line)
			line=str(' </ul>')+'\n'
			f.write(line)















			line=str("<head>")+'\n'
			f.write(line)
                        finder = SpacegroupAnalyzer(strt)
                        ##sgp=str(finder.get_spacegroup_number())
                        sgp=str(finder.get_spacegroup_symbol())
                        sgp.replace('/', 'slash')
                        print ("Checking ITssssssssssssssssssss 1L",file)
                        if '1L' in file:
                               typ=str('1L')
                               print ("ITssssssssssssssssssss 1L")
                        elif '2L' in file:
                               typ=str('2L')
                        elif 'Mol' in file:
                               typ=str('Mol')
                        else:
                               typ=str('Bulk')
			line=str('<a href="http://www.nist.gov/public_affairs/disclaimer.cfm"> NIST Disclaimer</a>')+'\n'
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Structural formula: ")+str(formula)+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Functional: ")+str(tag)+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Space group : ")+sgp+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Calculation type: ")+typ+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+ str(" JARVIS ID: ")+str(name)+str("</h2>")
			f.write(line)
                        if str(final_enp)!='na':
                           line=str('<h4 style="color:blue;">')+str("Final energy/atom (eV): ")+str(final_enp)+str("</h4>")
		    	   f.write(line)
                        line=str('<h4 style="color:blue;">')+ str(abc)+str("</h4>")
			f.write(line)
                        line=str('<h4 style="color:blue;">')+ str(angles)+str("</h4>")
			f.write(line)
                        path_to_file=str("VASP-FIGS/")+str(name)+str("/")+str("MAIN-RELAX.zip")
                        #line=str('<a href=')+str(path_to_file)+str(' >Download input files</a>')
                        line=str('<a href=')+str(path_to_file)+str(' download=')+str(name)+str('>Download input files</a>')
			f.write(line)
			line=str("</script>")+'\n'
			f.write(line)
			line=str('<script type="text/javascript"  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>')+'\n'
			f.write(line)
			line=str('<script type="text/javascript" src="JSmol.min.js"></script>')+'\n'
			f.write(line)
			line=str('<script type="text/javascript">')+'\n'
			f.write(line)
			line=str("Info = {")+'\n'
			f.write(line)
			line=str("	width: 400,")+'\n'
			f.write(line)
			line=str("	height: 400,")+'\n'
			f.write(line)
			line=str("	debug: false,")+'\n'
			f.write(line)
			line=str('	color: "0xC0C0C0",')+'\n'
			f.write(line)
			line=str("  disableJ2SLoadMonitor: false,")+'\n'
			f.write(line)
			line=str("  disableInitialConsole: true,")+'\n'
			f.write(line)
			line=str("	addSelectionOptions: false,")+'\n'
			f.write(line)
			line=str('	use: "HTML5",')+'\n'
			f.write(line)
			line=str("	readyFunction: null,")+'\n'
			f.write(line)
	#line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
			line=str('	script: "load VASP-FIGS/')+str(name)+str("/")+str(name)+str('_prem.cif"')+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("loadJmol = function() {")+'\n'
			f.write(line)
			line=str('  $("#mathdiv").show();')+'\n'
			f.write(line)
			line=str('	$("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))')+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("$(document).ready(function(){")+'\n'
			f.write(line)
			line=str("  checkMathJax(loadJmol);")+'\n'
			f.write(line)
			line=str("});")+'\n'
			f.write(line)
			line=str("checkMathJax = function(fReady) {")+'\n'
			f.write(line)
			line=str('  return (!$("#MathJax_Message")[0] ')+'\n'
			f.write(line)
			line=str("    || $(")+str("'")+str('#MathJax_Message:contains("Loading")')+str("'")+str(")[0]")+'\n'
			f.write(line)
			line=str("    ? setTimeout(function(){checkMathJax(fReady)},10)")+'\n'
			f.write(line)
			line=str("    : fReady());")+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("</script>")+'\n'
			f.write(line)
			line=str("</head>")+'\n'
			f.write(line)
			line=str("<body>")+'\n'
			f.write(line)
			line=str('<table><tr><td  style="width:450px">')+'\n'
			#line=str('<table><tr><td valign=top style="width:450px">')+'\n'
			#f.write(line)
			#line=str("</td><td>")+'\n'
			f.write(line)
			line=str("<span id=mydiv></span>")+'\n'
			f.write(line)
			line=str("</td></tr></table>")+'\n'
			f.write(line)



			#en_convg=EN_fig(pref='ENCUT',storedir=dir_name)
                        print "cwd",os.getcwd()
			#search,k_convg=K_fig(pref='KPOINTS',storedir=dir_name)
			#en_convg=EN_fig(pref='ENCUT',storedir=dir_name)
			try:
			   group=group+1
			   search,k_convg=K_fig(pref='KPOINTS',storedir=dir_name)

			   en_convg=EN_fig(pref='ENCUT',storedir=dir_name)
			   print ("search===",search,k_convg)
			   line=str(key)+str(',')+str(search)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           info.append(str(search))
                           info.append(str(name))
                           info.append(str(formula))
                           info.append(str(tag))
                           xc_info['search']=str(search)
                           xc_info['formula']=str(formula)
                           xc_info['name']=str(name)
                           xc_info['tag']=str(tag)
			   #print ("line======",line)
			   group=group+1
			   line=str(key)+str(',')+str(name)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           ref='NA'
                           typ='na'
                           print "GOING HERE 1"








			   ref_f=open("REFERENCE","r")
			   for rline in ref_f:
			       r=(str(rline))#.split('\n')[0]
                               if 'mp' in r and 'L' not in r:
                                  try:
                                     ref=(str(r)).split('\n')[0]
                                  except:
                                      pass
                               if '1L' in file: #str(l.split('.vasp')[0]).split('POSCAR-')[1]
                                   #r=str(rline.split('.vasp')[0]).split('POSCAR-')[1]
                                   #import sys
                                   #sys.exit()
                                   try:
                                      r=str(rline.split('.vasp')[0]).split('POSCAR-')[1]
			              #ref=str(str((str(line)).split('\n')[0]).split('.vasp')[0])split('POSCAR-')[1]
                                   except:
                                      pass
                               all_refs.append(str(r))
                           #icsd=get_struct_from_mp('mp-134')
                           icsd='na'
                           print "REFERENCE ISSSSSSSSS",all_refs,type(ref)
                           #icsd=get_struct_from_mp(ref)
                           try:
                              icsd=get_struct_from_mp(ref)
                           except:
                              pass
                           print "GOING HERE 2"
                           print "icsd=",icsd
			   group=group+1
			   line=str(key)+str(',')+str(ref)+str(',')+str(period)+str(',')+str(group)+'\n'
                           ################################################################################info.append(str(ref))
                           ##icsd='na'
                           if '1L' in file:
                               typ=str('1L')
                           elif '2L' in file:
                               typ=str('2L')
                           elif 'Mol' in file:
                               typ=str('Mol')
                           else:
                               typ=str('Bulk')
                           info.append(typ)
                           xc_info['typ']=str(typ)

			   #df.write(line)

			   line=str('<!--Source:')+str(os.getcwd())+str(' -->')+'\n'
			   f.write(line)
			   line=str('<h2 style="color:blue;">Convergence</h2>')+'\n'
			   f.write(line)
                           #line=str('<p > This is convergenece test </p>')+'\n'
                           line=str('<p > Calculations are done using VASP software. Convergence on KPOINTS and ENCUT is done with respect to total energy of the system within 0.001 eV tolerance. Please note convergence on KPOINTS and ENCUT is generally done for target properties, but here we assume energy-convergence with 0.001 eV should be sufficient for other properties also. The points on the curves are obtained with single-point calculation (nuber of ionic steps,NSW=1). However, for very accurate calculations, NSW>1 might be needed.   </p>')+'\n'
	                   f.write(line)
			   line=str('    <div class="row">')+'\n'
			   f.write(line)
			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Encut.png" width="450" height="350">')+'\n'
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
			   f.write(line)
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/Encut.png" width="450" height="350">')+'\n'
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('      </div>')+'\n'
			   f.write(line)
			except:
			    pass
			line=str('<h2 style="color:blue;">Structural analysis</h2>')+'\n'
			f.write(line)
                        line=str('<p > The following shows the X-ray diffraction (XRD) pattern and the Radial distribution function (RDF) plots. XRD peaks should be comparable to experiments. Relative intensities may differ.  </p>')+'\n'
                        #line=str('<p > This is Structural analysis </p>')+'\n'
	                f.write(line)
			line=str('    <div class="row">')+'\n'
			f.write(line)
			line=str('      <div class="col s12 m6">')+'\n'
			f.write(line)
			line=str('        <image src="VASP-FIGS/')+str(name)+str('/xrd.png" width="450" height="350">')+'\n'
			f.write(line)
			line=str('        <image src="VASP-FIGS/')+str(name)+str('/rdf.png" width="450" height="350">')+'\n'
			f.write(line)
			line=str('      </div>')+'\n'
			f.write(line)
                        bandgap=str('na')
			#bandgap=banddos(pref='MAIN-BAND',storedir=dir_name)
                        #print ("bandgappppppppppppppppppppppppppppppp",bandgap)
			#bandgap=banddos(pref='MAIN-BAND',storedir=dir_name)
			try:

			   bandgap=banddos(pref='MAIN-BAND',storedir=dir_name)
                           
                           print ("bandgappppppppppppppppppppppppppppppp",bandgap)
			except:
			    pass
			group=group+1
			line=str(key)+str(',')+str(bandgap)+str(',')+str(period)+str(',')+str(group)+'\n'
                        info.append(str(bandgap))
                        xc_info['bandgap']=str(bandgap)
			#df.write(line)



			line=str('<h2 style="color:blue;">Electronic structure</h2>')+'\n'
			f.write(line)
                        line=str('<p > The following shows the electronic density of states and bandstructure. DFT is generally predicted to underestimate bandgap of materials. Accurate band-gaps are obtained with higher level methods (with high computational requirement) such as HSE, GW, which are under progress. Total DOS, Orbital DOS and Element dos buttons are provided for density of states options. Energy is rescaled to make Fermi-energy zero. In the bandstructure plot, spin up is is shown with blue lines while spin down are shown with red lines. Non-degenerate spin-up and spin-down states (if applicable) would imply a net orbital magnetic moment in the system.  </p>')+'\n'
                        #line=str('<p > This is analysis </p>')+'\n'
	                f.write(line)
			line=str('<p>')+str('Bandgap (eV): ')+str(bandgap)+str('</p>')+'\n'
			f.write(line)
			line=str('      <div class="col s12 m6">')+'\n'
			f.write(line)
			line=str('        <image id="DOS" src="VASP-FIGS/')+str(name)+str('/TDos.png" width="450" height="350">')+'\n'
			f.write(line)
			line=str('        <image src="VASP-FIGS/')+str(name)+str('/Banddi.png" width="450" height="350">')+'\n'
			f.write(line)
			line=str('      </div>')+'\n'
			f.write(line)
			line=str('<button onclick="myTDOSFunction()">Total DOS</button>')+'\n'
			f.write(line)

			line=str('<script>')+'\n'
			f.write(line)
			line=str('function myTDOSFunction() {')+'\n'
			f.write(line)
			line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/TDos.png" ;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)
			line=str('</script>')+'\n'
			f.write(line)
			line=str('<button onclick="myPDOSFunction()">Orbital DOS</button>')+'\n'
			f.write(line)

			line=str('<script>')+'\n'
			f.write(line)
			line=str('function myPDOSFunction() {')+'\n'
			f.write(line)
			line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/Dos.png" ;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)
			line=str('</script>')+'\n'
			f.write(line)
			line=str('<button onclick="myEDOSFunction()">Element DOS</button>')+'\n'
			f.write(line)

			line=str('<script>')+'\n'
			f.write(line)
			line=str('function myEDOSFunction() {')+'\n'
			f.write(line)
			line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/EDos.png" ;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)
			line=str('</script>')+'\n'
			f.write(line)



			#bandgapHSE=banddosHSE(pref='MAIN-HSE',storedir=dir_name)
                        #print ("HSE DONE",dir_name,bandgapHSE)
                        #import sys
                        #sys.exit()
                        ccwd1=str(os.getcwd())
                        dif=str('na')
                        try:
			                dif=vtot(pref='MAIN-RELAX',storedir=storedir,name=str(name))
					line=str('<h2 style="color:blue;">Electrostatic potential</h2>')+'\n'
					f.write(line)
					line=str('<p > The following plot shows the plane averaged electrostatic potential (ionic+Hartree) along z direction. The red line shows the Fermi-energy while the green line shows the maximum value of the electrostatic potential. For slab structures (with vacuum along z-direction), the difference in these two values can be used to calculate work-function of the material.   </p>')+'\n'
					#line=str('<p > This is analysis </p>')+'\n'
					f.write(line)
					line=str('      <div class="col s12 m6">')+'\n'
					f.write(line)
					line=str('        <image id="DOS" src="VASP-FIGS/')+str(name)+str('/ESTAT.png" width="450" height="350">')+'\n'
					f.write(line)
					line=str('      </div>')+'\n'
					f.write(line)
                        except:
                           pass
                        os.chdir(ccwd1)
                        xc_info['dif']=str(dif)





                        bandgapHSE='na'
			try:
                           dir_list=[d for d in os.listdir('./') if os.path.isdir(os.path.join('./', d))]
                           for dl in dir_list:
                               if "MAIN-HSE" in dl:
                                   try:
                                     HSEdone=(Vasprun(str(os.getcwd())+str("/")+str(dl)+str("/")+str("vasprun.xml"))).converged
                                   except:
                                      pass
                                   if HSEdone==True:
                                        print "HSE folder=True",str(os.getcwd())+str("/")+str(dl)+str("/")+str("vasprun.xml")
			                bandgapHSE=banddosHSE(pref='MAIN-HSE',storedir=dir_name)
                           
                                        print ("bandgapppppppppppppppppppppppppppppppHSE",bandgapHSE)
					group=group+1
					#line=str(key)+str(',')+str(bandgap)+str(',')+str(period)+str(',')+str(group)+'\n'
					#info.append(str(bandgap))
					#df.write(line)



					line=str('<h2 style="color:blue;">HSE Electronic structure</h2>')+'\n'
					f.write(line)
					line=str('<p > The following shows the electronic density of states and bandstructure using HSE06 hybrid functional.  </p>')+'\n'
					#line=str('<p > This is analysis </p>')+'\n'
					f.write(line)
			                line=str('<p>')+str('Bandgap (eV): ')+str(bandgapHSE)+str('</p>')+'\n'
			                f.write(line)
					line=str('      <div class="col s12 m6">')+'\n'
					f.write(line)
					line=str('        <image id="DOS" src="VASP-FIGS/')+str(name)+str('/TDosHSE.png" width="450" height="350">')+'\n'
					f.write(line)
					line=str('        <image src="VASP-FIGS/')+str(name)+str('/BanddiHSE.png" width="450" height="350">')+'\n'
					f.write(line)
					line=str('      </div>')+'\n'
					f.write(line)
					line=str('<button onclick="myTDOSFunction()">Total DOS</button>')+'\n'
					f.write(line)

					line=str('<script>')+'\n'
					f.write(line)
					line=str('function myTDOSFunction() {')+'\n'
					f.write(line)
					line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/TDosHSE.png" ;')+'\n'
					f.write(line)
					line=str('}')+'\n'
					f.write(line)
					line=str('</script>')+'\n'
					f.write(line)
					line=str('<button onclick="myPDOSFunction()">Orbital DOS</button>')+'\n'
					f.write(line)

					line=str('<script>')+'\n'
					f.write(line)
					line=str('function myPDOSFunction() {')+'\n'
					f.write(line)
					line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/DosHSE.png" ;')+'\n'
					f.write(line)
					line=str('}')+'\n'
					f.write(line)
					line=str('</script>')+'\n'
					f.write(line)
					line=str('<button onclick="myEDOSFunction()">Element DOS</button>')+'\n'
					f.write(line)

					line=str('<script>')+'\n'
					f.write(line)
					line=str('function myEDOSFunction() {')+'\n'
					f.write(line)
					line=str('    document.getElementById("DOS").src = "VASP-FIGS/')+str(name)+str('/EDosHSE.png" ;')+'\n'
					f.write(line)
					line=str('}')+'\n'
					f.write(line)
					line=str('</script>')+'\n'
					f.write(line)





			except:
			    pass
                        info.append(str(bandgapHSE))
                        xc_info['bandgapHSE']=str(bandgapHSE)










			#optics(pref='MAIN-OPTICS',storedir=dir_name)
			try:
			   line=str('<h2 style="color:blue;">Optical properties </h2>')+'\n'
			   f.write(line)
                           line=str('<p > Incident photon energy dependence of optical is shown below. Only interband optical transitions are taken into account.Please note the underestimatation of band-gap problem with DFT will reflect in the spectra as well. For very accurate optical properties GW/BSE calculation would be needed, which is yet to be done because of their very high computational cost. Optical properties for layered materials needs to be rescaled with the actual thickness to simulation z-box ratio. Absorption coeffiecient is in cm-1 unit. </p>')+'\n'
                           #line=str('<p > This is analysis </p>')+'\n'
	                   f.write(line)
			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   optics(pref='MAIN-OPTICS',storedir=dir_name)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Real.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Imag.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str("      </div>")+'\n'
			   f.write(line)


			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Absorp.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/ELS.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str("      </div>")+'\n'
			   f.write(line)



			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Refr.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Extc.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str("      </div>")+'\n'
			   f.write(line)


			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Opt_con" width="450" height="350">')+'\n'
			   f.write(line)
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/Extc.png" width="450" height="350">')+'\n'
			   #f.write(line)
			   line=str("      </div>")+'\n'
			   f.write(line)

			except:
			    pass
			#c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,ph_modes,pg_phonopy=elastic_tensor(pref='MAIN-ELASTIC',storedir=dir_name,name=str("JV2D-")+str(count))
			#BV=float((c11+c22+c33)+2*(c12+c23+c31))/float(9)
			#GV=float((c11+c22+c33)-(c12+c23+c31)+3*(c44+c55+c66))/float(15)
                        BV='na'
                        GV='na'
			#c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,ph_modes,pg_phonopy=elastic_tensor(pref='MAIN-ELASTIC',storedir=dir_name,name=str("JVASP-")+str(count))
			#c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,ph_modes,pg_phonopy=elastic_tensor(pref='MAIN-ELASTIC',storedir=dir_name,name=str("JVASP-")+str(count))
			try:
			   c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,ph_modes,pg_phonopy=elastic_tensor(pref='MAIN-ELASTIC',storedir=dir_name,name=name)
			   line=str('<h2 style="color:blue;">Elastic tensor and derived phonon properties</h2>')+'\n'
			   f.write(line)
                           line=str('<p > Elastic tensor calculated for the conventional cell of the system with finite-difference method. For layered materials, the elastic constants are rescaled with respect to vacuum padding (see the input files) and the units for elastic coefficients are in N/m. Phonons obtained from this calcuation are also shown.</p>')+'\n'
	                   f.write(line)
                           line=str('<p > WARNING: Please note this may not be the exact phonon modes of the system as we did not test the cell-size dependence of phonons yet. At least 1.2 nm x1.2 nm x1.2 nm or more is needed for obtaining reliable phonon spectrum. For systems having primitive-cell phonon representation tables, I denotes infrared activity and R denotes Raman active modes (where applicabale).  </p>')+'\n'
                           #line=str('<p > This is analysis </p>')+'\n'
	                   f.write(line)
			   BV=float((c11+c22+c33)+2*(c12+c23+c31))/float(9)
			   GV=float((c11+c22+c33)-(c12+c23+c31)+3*(c44+c55+c66))/float(15)
			   BV=round(BV,3)
			   GV=round(GV,3)
			   group=group+1
			   line=str(key)+str(',')+str(BV)+str(',')+str(period)+str(',')+str(group)+'\n'
                           #info.append(str(ref))
			   #df.write(line)

			   group=group+1
			   line=str(key)+str(',')+str(GV)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           line=str("<h4>Elastic Tensor C<sub>ij</sub> GPa</h4> ")#+'\n'
			   f.write(line)

			   line=str('    <div class="row">')+'\n'
			   f.write(line)
			   line=str('      <div class="span5">')+'\n'
			   f.write(line)
			   line=str('        <div class="matrix">')+'\n'
			   f.write(line)
			   line=str('          <table  class="table">')+'\n'
			   f.write(line)
			   line=str("            <thead>")+'\n'
			   f.write(line)
			   line=str('              <th colspan="6" >')+'\n'
			   f.write(line)
			   line=str('                <h4 class="text-center">')+'\n'
			   f.write(line)
			   #line=str('Elastic Tensor C<sub>ij</sub> GPa')+'\n'
			   #f.write(line)
			   #line=str('<span class="units"">(GPa)</span>')+'\n'
			   #f.write(line)
			   line=str('                </h4>')+'\n'
			   f.write(line)
			   line=str("              </th>")+'\n'
			   f.write(line)
			   line=str("            </thead>")+'\n'
			   f.write(line)
			   line=str("            <tbody>")+'\n'
			   f.write(line)

                           #line=str("<h4>Bulk Modulus B<sub>V</sub></h4>")+'\n'
			   #f.write(line)
                           #line=str("<p> ")+str(Bv)+str(" GPa </p>")+'\n'
			   #f.write(line)
                           #line=str("<h4>Shear Modulus G<sub>V</sub></h4>")+'\n'
			   #f.write(line)
                           #line=str("<p> ")+str(Gv)+str(" GPa </p>")+'\n'
			   #f.write(line)
			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c11)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c12)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c13)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c14)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c15)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c16)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c21)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c22)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c23)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c24)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c25)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c26)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c31)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c32)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c33)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c34)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c35)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c36)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c41)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c42)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c43)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c44)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c45)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c46)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)



			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c51)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c52)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c53)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c54)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c55)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c56)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c61)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c62)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c63)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c64)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c65)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(c66)+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)

			   line=str("            </tbody>")+'\n'
			   f.write(line)
			   line=str("          </table>")+'\n'
			   f.write(line)
			   line=str("        </div>")+'\n'
			   f.write(line)
			   line=str("      </div>")+'\n'
			   f.write(line)
			   line=str("<table>")+'\n'
			   f.write(line)
			   line=str("  <tr>")+'\n'
			   f.write(line)
			   line=str("    <th>Phonon mode (cm<sup>-1</sup>)</th>")+'\n'
			   f.write(line)
			   for modes in ph_modes:
				   line=str("                <tr>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(modes)+str("</td>")+'\n'
				   f.write(line)
				   line=str("                </tr>")+'\n'
				   f.write(line)
			   line=str("</table>")+'\n'
			   f.write(line)
			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/pdos_conv.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/PBAND.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('      </div>')+'\n'
			   f.write(line)
			   line=str("    </div>")+'\n'
			   f.write(line)
                           xc_info['elastic_prop']=list([ c11,c12,c13,c14,c15,c16,c21,c22,c23,c24,c25,c26,c31,c32,c33,c34,c35,c36,c41,c42,c43,c44,c45,c46,c51,c52,c53,c54,c55,c56,c61,c62,c63,c64,c65,c66,list(ph_modes),pg_phonopy])
			except:
			      pass
                        #line=str("<h4>Bulk Modulus B<sub>V</sub></h4>")#+'\n'
			#f.write(line)
                        #line=str("<p> ")+str(BV)+str(" GPa </p>")#+'\n'
			#f.write(line)
                        #line=str("<h4>Shear Modulus G<sub>V</sub></h4>")#+'\n'
			#f.write(line)
                        #line=str("<p> ")+str(GV)+str(" GPa </p>")#+'\n'
			#f.write(line)
                        info.append(str(BV))
                        info.append(str(GV))
                        xc_info['BV']=str(BV)
                        xc_info['GV']=str(GV)
                         
                        cwd1=str(os.getcwd())
			#nmodes=phonon(storedir=storedir)
                        #print ("nmodes=",nmodes)
			try:
			   nmodes=phonon(storedir=storedir,name=name)


			   #line=str("<h2>Phonon</h2>")+'\n'
			   #f.write(line)
			   line=str("<h4>Point group</h4>")+'\n'
			   f.write(line)
			   line=str("<p>")+str(pg_phonopy)+str("</p>")+'\n'
			   f.write(line)
			   line=str("<table>")+'\n'
			   f.write(line)
			   line=str("  <tr>")+'\n'
			   f.write(line)
			   line=str("    <th>Phonon mode (cm<sup>-1</sup>)</th>")+'\n'
			   f.write(line)
			   #line=str(" Visualize Phonos")+'\n'
                           #purl=http://www.ctcms.nist.gov/~knc6/jsmol/VASP-FIGS/JVASP-60/EDos1.png
                           #purl=str("VASP-FIGS/")+str("JVASP-")+str(count)+str("/")+str("vphonon.html")
    #vfile=str(storedir)+str("/")+str(name)+str("_vphonon.html")
                           #purl=str("http://www.ctcms.nist.gov/~knc6/jsmol/")+str("JVASP-")+str(count)+str("/")+str("vphonon.html")
                           #purl=str("JVASP-")+str(count)+str("_vphonon.html")
                           #purl=str("JVASP-")+str(name)+str("_vphonon.html")
                           purl=str(name)+str("_vphonon.html")
                           line=str('<a href=')+str(purl)+str(' target="_blank" >Visualize Phonons here</a>')
			   f.write(line)
			   line=str("    <th>Representation</th>")+'\n'
			   f.write(line)
			   for modes in nmodes:
			       try:
				   line=str("                <tr>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(float(modes['frequency']))+str("</td>")+'\n'
				   #line=str("                    <td>")+str(float(33.35641)*float(modes['frequency']))+str("</td>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(modes['ir_label'])+str("</td>")+'\n'
				   f.write(line)
				   line=str("                </tr>")+'\n'
				   f.write(line)
			       except:
				   line=str("                <tr>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(modes['frequency'])+str("</td>")+'\n'
				   f.write(line)
				   line=str("                </tr>")+'\n'
				   f.write(line)
			   line=str("  </tr>")+'\n'
			   f.write(line)
			   line=str("</table>")+'\n'
			   f.write(line)
			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/pdos_conv.png" width="450" height="350">')+'\n'
			   #f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/PBAND.png" width="450" height="350">')+'\n'
			   #f.write(line)
			   line=str('      </div>')+'\n'
			   f.write(line)
			   line=str("    </div>")+'\n'
			   f.write(line)


			except:
			   pass
                        #os.chdir(cwd1)
                        #cwd1=str(os.getcwd())
                        fin_en='na'
                        print ("cwd=",os.getcwd())
			#fin_en,search,num,avg,additional_info=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
                        os.chdir(ccwd1)
                        import sys
                        print ("sys")
			#fin_en,search,num,avg,additional_info=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
			#fin_en,search,num,avg,additional_info=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=name)
			try:
			   fin_en,search,num,avg,additional_info=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=name)
                           os.chdir(ccwd1)

                        #print "Effective mass",avg['p'][300][0]
                        #print "Effective mass",avg['p'][300][0][0][0],avg['p'][300][0][0][1],avg['p'][300][0][0][2]
                        #print "Effective mass",avg['p'][300][0][1][0],avg['p'][300][0][1][1],avg['p'][300][0][1][2]
                        #print "Effective mass",avg['p'][300][0][2][0],avg['p'][300][0][2][1],avg['p'][300][0][2][2]
			   group=group+1
			   line=str(key)+str(',')+str(fin_en)+str(',')+str(period)+str(',')+str(group)+'\n'
			   group=group+1
			   line=str(key)+str(',')+str(tag)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           info.append(str(num))
                           xc_info['SG']=str(num)
                           xc_info['CWD']=str(os.getcwd())



			   line=str('<h2 style="color:blue;">Thermoelectric properties</h2>')+'\n'

                           print ("Effective mass",avg['n'][300][0])
			   #line=str('<h2 style="color:blue;">Thermoelectric properties</h2>')+'\n'
			   f.write(line)
                           #line=str('<p > This is analysis </p>')+'\n'
                           line=str('<p > Thermoelectric properties are calculated using BoltzTrap code. Electron and hole mass tensors are given at 300 K. Following plots show the Seebeck coefficient and ZT factor (eigenvalues of the tensor shown) at 300 K along three different crystallographic directions. Seebeck coefficient and ZT plots can be compared for three different temperatures available through the buttons given below. Generally very high Kpoints are needed for obtaining thermoelectric properties. We assume the Kpoints obtained from above convergence were sufficient. </p>')+'\n'
	                   f.write(line)
                           line=str('<p > WARNING: Constant relaxation time approximation (10<sup>-14</sup> s) and only electronic contribution to thermal conductivity  were utilized for calculating ZT. </p>')+'\n'
	                   f.write(line)
                           line=str("<h3>Electron mass tensor (m<sub>e</sub> unit)</h3>")#+'\n'
			   f.write(line)
			   line=str('    <div class="row">')+'\n'
			   f.write(line)
			   line=str('      <div class="span5">')+'\n'
			   f.write(line)
			   line=str('        <div class="matrix">')+'\n'
			   f.write(line)
			   line=str('          <table  class="table">')+'\n'
			   f.write(line)
			   line=str("            <thead>")+'\n'
			   f.write(line)
			   line=str('              <th colspan="6" >')+'\n'
			   f.write(line)
			   line=str('<h3 >')+'\n'
			   #line=str('<h3 style="color:blue;" >')+'\n'
			   #f.write(line)
			   #line=str('Electron mass Tensor (m<sup>e</sup> unit) ')+'\n'
			   #f.write(line)
			   line=str('</h3>')+'\n'
			   f.write(line)
			   line=str("              </th>")+'\n'
			   f.write(line)
			   line=str("            </thead>")+'\n'
			   f.write(line)
			   line=str("            <tbody>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][0][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][0][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][0][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)

			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][1][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][1][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][1][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][2][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][2][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['n'][300][0][2][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)
			   line=str("            </tbody>")+'\n'
			   f.write(line)
			   line=str("          </table>")+'\n'
			   f.write(line)
			   line=str("        </div>")+'\n'
			   f.write(line)



                           line=str("<h3>Hole mass tensor  (m<sub>e</sub> unit)</h3>")#+'\n'
			   f.write(line)
			   line=str('    <div class="row">')+'\n'
			   f.write(line)
			   line=str('      <div class="span5">')+'\n'
			   f.write(line)
			   line=str('        <div class="matrix">')+'\n'
			   f.write(line)
			   line=str('          <table  class="table">')+'\n'
			   f.write(line)
			   line=str("            <thead>")+'\n'
			   f.write(line)
			   line=str('              <th colspan="6" >')+'\n'
			   f.write(line)
			   line=str('<h3 >')+'\n'
			   #line=str('<h3 style="color:blue;" >')+'\n'
			   f.write(line)
			   #line=str('Hole mass Tensor ')+'\n'
			   #f.write(line)
			   line=str('</h3>')+'\n'
			   f.write(line)
			   line=str("              </th>")+'\n'
			   f.write(line)
			   line=str("            </thead>")+'\n'
			   f.write(line)
			   line=str("            <tbody>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][0][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][0][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][0][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)

			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][1][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][1][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][1][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)


			   line=str("                <tr>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][2][0])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][2][1])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                    <td>")+str(avg['p'][300][0][2][2])+str("</td>")+'\n'
			   f.write(line)
			   line=str("                </tr>")+'\n'
			   f.write(line)
			   line=str("            </tbody>")+'\n'
			   f.write(line)
			   line=str("          </table>")+'\n'
			   f.write(line)
			   line=str("        </div>")+'\n'

                           xc_info.update(avg)#list(avg['n'][300])
                           #xc_info['p']=list(avg['p'][300])

			   line=str('<div class="row">')+'\n'
			   f.write(line)
			   line=str('    <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('      <image id="seebeck" src="VASP-FIGS/')+str(name)+str('/S_mu_300.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('      <image id="zt" src="VASP-FIGS/')+str(name)+str('/ZT_mu_300.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('    </div>')+'\n'
			   f.write(line)




			   line=str('<button onclick="myFunction1()">Temp-100 K</button>')+'\n'
			   f.write(line)

			   line=str('<script>')+'\n'
			   f.write(line)
			   line=str('function myFunction1() {')+'\n'
			   f.write(line)
			   line= str('  document.getElementById("seebeck").src = "VASP-FIGS/')+str(name)+str('/S_mu_100.png" ;')+'\n'
			   f.write(line)
			   line=str('  document.getElementById("zt").src ="VASP-FIGS/')+str(name)+str('/ZT_mu_100.png";')+'\n'
			   f.write(line)
			   line=str('}')+'\n'
			   f.write(line)
			   line=str('</script>')+'\n'
			   f.write(line)

			   line=str('<button onclick="myFunction2()">Temp-300 K</button>')+'\n'
			   f.write(line)

			   line=str('<script>')+'\n'
			   f.write(line)
			   line=str('function myFunction2() {')+'\n'
			   f.write(line)
			   line=str('    document.getElementById("seebeck").src = "VASP-FIGS/')+str(name)+str('/S_mu_300.png" ;')+'\n'
			   f.write(line)
			   line=str('    document.getElementById("zt").src = "VASP-FIGS/')+str(name)+str('/ZT_mu_300.png";')+'\n'
			   f.write(line)
			   line=str('}')+'\n'
			   f.write(line)
			   line=str('</script>')+'\n'
			   f.write(line)


			   line=str('<button onclick="myFunction3()">Temp-600 K</button>')+'\n'
			   f.write(line)

			   line=str('<script>')+'\n'
			   f.write(line)
			   line=str('function myFunction3() {')+'\n'
			   f.write(line)
			   line=str('    document.getElementById("seebeck").src =  "VASP-FIGS/')+str(name)+str('/S_mu_600.png";')+'\n'
			   f.write(line)
			   line=str('    document.getElementById("zt").src = "VASP-FIGS/')+str(name)+str('/ZT_mu_600.png";')+'\n'
			   f.write(line)
			   line=str('}')+'\n'
			   f.write(line)
			   line=str('</script>')+'\n'
			   f.write(line)



			   #line=str('<button onclick="myFunction3a()">Rel-10^-12</button>')+'\n'
			   #f.write(line)

			   #line=str('<script>')+'\n'
			   #f.write(line)
			   #line=str('function myFunction3a() {')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("seebeck").src =  "VASP-FIGS/')+str(name)+str('/S_mu_300_rel_12.png";')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("zt").src = "VASP-FIGS/')+str(name)+str('/ZT_mu_300_rel_12.png";')+'\n'
			   #f.write(line)
			   #line=str('}')+'\n'
			   #f.write(line)
			   #line=str('</script>')+'\n'
			   #f.write(line)



			   #line=str('<button onclick="myFunction3b()">Rel-10^-14</button>')+'\n'
			   #f.write(line)

			   #line=str('<script>')+'\n'
			   #f.write(line)
			   #line=str('function myFunction3b() {')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("seebeck").src =  "VASP-FIGS/')+str(name)+str('/S_mu_300.png";')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("zt").src = "VASP-FIGS/')+str(name)+str('/ZT_mu_300.png";')+'\n'
			   #f.write(line)
			   #line=str('}')+'\n'
			   #f.write(line)
			   #line=str('</script>')+'\n'
			   #f.write(line)


			   #line=str('<button onclick="myFunction3c()">Rel-10^-16</button>')+'\n'
			   #f.write(line)

			   #line=str('<script>')+'\n'
			   #f.write(line)
			   #line=str('function myFunction3c() {')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("seebeck").src =  "VASP-FIGS/')+str(name)+str('/S_mu_300_rel_16.png";')+'\n'
			   #f.write(line)
			   #line=str('    document.getElementById("zt").src = "VASP-FIGS/')+str(name)+str('/ZT_mu_300_rel_16.png";')+'\n'
			   #f.write(line)
			   #line=str('}')+'\n'
			   #f.write(line)
			   #line=str('</script>')+'\n'
			   #f.write(line)









			   line=str('<h2 style="color:blue;">Reference</h2>')+'\n'
			   f.write(line)
                           for r in all_refs:
                               if 'mp' in r and 'L' not in r:
                                  try:
                                     ref=(str(r)).split('\n')[0]
                                  except:
                                      pass
                               line=str('<br>')+str(r)+'<br>'
			       f.write(line)
                           
                           if icsd !='na':
                        
                             line=str('<p>ICSD-ID: ')+str(icsd)+str('</p>')+'\n'
			     f.write(line)
                             comp=strt.composition#red_formula.split(' ')
                             factor = gcd(*(int(i) for i in comp.values()))
                             comp1=comp/factor
                             tmp=comp1.alphabetical_formula
                             new_formula=''.join(tmp).replace(" ","")
                             url=str('http://aflowlib.org/material.php?id=')+str(new_formula)+str("_")+str("ICSD")+str("_")+str(icsd)
                             try:
                               r = requests.get(url)
                               web_content=r.text.splitlines()
                               aflow=False
                               for z in web_content:
                                  if 'ICSD' in z:
                                       aflow=True
                               if aflow==True and icsd!=None and icsd !=[]:
                                  
                                 line=str('<a href=')+str(url)+str(' target="_blank" >AFLOW link</a>')
			         f.write(line)
                                 line=str('<br>')+'\n'
			         f.write(line)
                                  
                             except: 
                                pass
                           print "ICSD ID IS=",icsd
                           if 'mp' in ref and 'L' not in ref:
                             link=str('https://www.materialsproject.org/materials/')+str(ref)
                             line=str('<br>')+str('<a href=')+str(link)+str(' target="_blank" >MP link</a>')
			     f.write(line)
                             line=str('<br>')+'\n'
			     f.write(line)

			   line=str("<!--")+'\n'
			   f.write(line)
			   line=str("  <tr>")+'\n'
			   f.write(line)
			   line=str("    <th>Parameters</th>")+'\n'
			   f.write(line)
			   line=str("    <th>Value</th>")+'\n'
			   f.write(line)
                           #xc_info['JARVIS-ID']=str(name)
                           vars=['PSCENC','TEWEN', 'DENC','EXHF','XCENC','PAW_DB','EENTRO','EBANDS','EATOM','mag_mom','NELECT']
                           #xc_info['INFO']=vars
                           #vars=['PSCENC','TEWEN', 'DENC','EXHF','XCENC','PAW_DB','EENTRO','EBANDS','EATOM','TOTEN','mag_mom','NELECT']
			   try:
			      for a,b in zip(vars,additional_info):
                                   xc_info[a]=b
                                   #print "a and b are", a,b
				   line=str("                <tr>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(a)+str("</td>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(b)+str("</td>")+'\n'
				   f.write(line)
				   line=str("                </tr>")+'\n'
				   f.write(line)
                              ##xc_info_arr.append(xc_info)
			   except:
                                    pass
			   line=str("  </tr>")+'\n'
			   f.write(line)
			   line=str("</table>")+'\n'
			   f.write(line)

                           line=str('<p>')+str(ref)+str('</p>')+'\n'
			   f.write(line)

			   line=str("-->")+'\n'
			   f.write(line)

			   #line=str('      <div class="col s12 m6">')+'\n'
			   #f.write(line)
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/S_mu_100.png" width="450" height="350">')+'\n'
			   #f.write(line)
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/ZT_mu_100.png" width="450" height="350">')+'\n'
			   #f.write(line)
			   #line=str('      </div>')+'\n'
			   #f.write(line)
			   #line=str("    </div>")+'\n'
			   #f.write(line)
			   line=str("</body>")+'\n'
			   f.write(line)
			   line=str("</html>")+'\n'
			   f.write(line)
			   f.close()
			except:
			   line=str('<h2 style="color:blue;">Reference</h2>')+'\n'
			   f.write(line)
                           for r in all_refs:
                               if 'mp' in r and 'L' not in r:
                                  try:
                                     ref=(str(r)).split('\n')[0]
                                  except:
                                      pass
                               line=str('<br>')+str(r)+'<br>'
			       f.write(line)
                           ##for r in all_refs:
                           ##       line=str(r)+'<br>'
			   ##       f.write(line)
                           #line=str(ref)+'\n'
			   #f.write(line)

                           if icsd !='na':
                        
                             line=str('<p>ICSD-ID: ')+str(icsd)+str('</p>')+'\n'
			     f.write(line)
                             comp=strt.composition#red_formula.split(' ')
                             factor = gcd(*(int(i) for i in comp.values()))
                             comp1=comp/factor
                             tmp=comp1.alphabetical_formula
                             new_formula=''.join(tmp).replace(" ","")
                             url=str('http://aflowlib.org/material.php?id=')+str(new_formula)+str("_")+str("ICSD")+str("_")+str(icsd)
                             try:
                               r = requests.get(url)
                               web_content=r.text.splitlines()
                               aflow=False
                               for z in web_content:
                                  if 'ICSD' in z:
                                       aflow=True
                               #if aflow==True:
                               if aflow==True and icsd!=None:
                                 line=str('<a href=')+str(url)+str(' target="_blank" >AFLOW link</a>')
			         f.write(line)
                                 line=str('<br>')+'\n'
			         f.write(line)
                                  
                             except: 
                                pass
                           if 'mp' in ref and 'L' not in ref:
                             link=str('https://www.materialsproject.org/materials/')+str(ref)
                             line=str('<a href=')+str(link)+str(' target="_blank" >MP link</a>')
			     f.write(line)
                             line=str('<br>')+'\n'
			     f.write(line)










                           
                     



			   line=str("<!--")+'\n'
			   f.write(line)
			   line=str("<table>")+'\n'
			   f.write(line)
			   line=str("  <tr>")+'\n'
			   f.write(line)
			   line=str("    <th>Parameters</th>")+'\n'
			   f.write(line)
			   line=str("    <th>Value</th>")+'\n'
			   f.write(line)
                           vars=['PSCENC','TEWEN', 'DENC','EXHF','XCENC','PAW_DB','EENTRO','EBANDS','EATOM','mag_mom','NELECT']
                           #vars=['PSCENC','TEWEN', 'DENC','EXHF','XCENC','PAW_DB','EENTRO','EBANDS','EATOM','TOTEN','mag_mom','NELECT']
			   try:
			      for a,b in zip(vars,additional_info):
                                   #print "a and b are", a,b
				   line=str("                <tr>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(a)+str("</td>")+'\n'
				   f.write(line)
				   line=str("                    <td>")+str(b)+str("</td>")+'\n'
				   f.write(line)
				   line=str("                </tr>")+'\n'
				   f.write(line)
			   except:
                                    pass
			   line=str("-->")+'\n'
			   f.write(line)
			   line=str("  </tr>")+'\n'
			   f.write(line)
			   line=str("</table>")+'\n'
			   f.write(line)
			   line=str("</body>")+'\n'
			   f.write(line)
			#line=str("</body>")+'\n'
			   line=str("</html>")+'\n'
			   f.write(line)
			   f.close()
			   pass
                        try:
                            info.append(str(float(fin_en)/float(len(strt))))
                            xc_info['enp']=str(str(float(fin_en)/float(len(strt))))
                        except:
                            xc_info['enp']=str('na')
                            info.append('na')
			#write_html(storedir=dir_name,name=str("JV-")+str(count),nmodes=nmodes,c11=c11,c12=c12,c13=c13,c14=c14,c15=c15,c16=c16,c21=c21,c22=c22,c23=c23,c24=c24,c25=c25,c26=c26,c31=c31,c32=c32,c33=c33,c34=c34,c35=c35,c36=c36,c41=c41,c42=c42,c43=c43,c44=c44,c45=c45,c46=c46,c51=c51,c52=c52,c53=c53,c54=c54,c55=c55,c56=c56,c61=c61,c62=c62,c63=c63,c64=c64,c65=c65,c66=c66)

			print(" ##MESSAGE##")
			print ("COPY JVASP-* folders in /u/WWW/knc6/jsmol/VASP-FIGS")
			print ("COPY JVASP-*.html folders in /u/WWW/knc6/jsmol")
			#print ("COPY JVASP-*_prem.cif folders in /u/WWW/knc6/jsmol/data")
        if done==False and partial_done==True :
                        formula=p_strt.composition.formula
			new_cal.append(file.split('_')[0])
			for cc in new_cal:
			   ccc=0
			   if cal in new_cal:
			      ccc=ccc+1

		  # for file in glob.glob('*_mp*'):
			group=0
			period=ccc

			print (file)
			count=count+1
			key=key+1
			#storedir=str('/data/knc6/NIST2/OLD/STORE')
			#storedir=str('/data/knc6/STORAGE/NEW_STORAGE')
			storedir=str('/rk2/knc6/FRESH32g_DIR')
                        with open('JARVIS-ID','r') as jf:
                             jl=jf.read()
                             name=jl.split('\n')[0] #str("JVASP-")+str(count)
			dir_name=os.path.join(storedir,str(name))
                        if not os.path.exists(dir_name): os.makedirs(dir_name)
			#fin_en,search,num,avg=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
			#fin_en,search,num,avg=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JVASP-")+str(count))
                        #print "Effective mass",avg['p'][300][0]
                        #print "Effective mass",avg['p'][300][0][0][0],avg['p'][300][0][0][1],avg['p'][300][0][0][2]
                        #print "Effective mass",avg['p'][300][0][1][0],avg['p'][300][0][1][1],avg['p'][300][0][1][2]
                        #print "Effective mass",avg['p'][300][0][2][0],avg['p'][300][0][2][1],avg['p'][300][0][2][2]
			print "COUNTTTT===",count,dir_name
                        ref_numb='na'
			#fin_en,search=boltz_bader(pref='MAIN-RELAX',storedir=storedir,name=str("JV2D-")+str(count))
			#if  os.path.exists(cwd+'/'+str(file)+str("/")+str("REFERENCE")): 
                        #   ref=open(cwd+'/'+str(file)+str("/")+str("TAGS")+str("REFERENCE"),"r")
                        #   for line in ref:
                        #       ref_numb=(str(line)).split('\n')[0]
                        #   ref.close()
			   #shutil.copy2(cwd+'/'+str(file)+str("/")+str("REFERENCE"),dir_name)
			tag='na'
                        print "tag",cwd+str('/')+str(file)+str("/")+str("FUNCTIONAL")
			if  os.path.exists("FUNCTIONAL"):
			#if  os.path.exists(cwd+str('/')+str(file)+str("/")+str("FUNCTIONAL")):
			    tag_file= open("FUNCTIONAL","r")
			    for line in tag_file:
				tag=(str(line)).split('\n')[0]
                                print ("tag====",tag)
			    tag_file.close()
			if tag !='na':
			   print ("New tag assigned")
			else:
			  tag=str("Unknown")
                        print "TAGGGG IS ",tag
			#name=str("JVASP-")+str(count)
                        filename=str(storedir)+str("/")+str(name)+str("/")+str(name)+str('_prem.cif')
                        p_strt.to(fmt= "cif", filename= filename)
			f=open(str(storedir)+str("/")+str(name)+str(".html"),"w")
			line=str("<!DOCTYPE html>")+'\n'
			f.write(line)
			line=str("<html>")+'\n'
			f.write(line)
			line=str('<body style="background-color:lightgreen;">')+'\n'
			f.write(line)
			line=str("<title>JARVIS-DFT</title>")+'\n'
			#line=str("<title>JSmol-MathJax Compatibility Test</title>")+'\n'
			f.write(line)




			line=str('<style>')+'\n'
			f.write(line)
			line=str('ul {')+'\n'
			f.write(line)
			line=str('    list-style-type: none;')+'\n'
			f.write(line)
			line=str('    margin: 0;')+'\n'
			f.write(line)
			line=str('    padding: 0;')+'\n'
			f.write(line)
			line=str('    overflow: hidden;')+'\n'
			f.write(line)
			line=str('    background-color: 	#808080;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)

			line=str('li {')+'\n'
			f.write(line)
			line=str('    float: left;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)

			line=str('li a {')+'\n'
			f.write(line)
			line=str('    display: inline-block;')+'\n'
			f.write(line)
			line=str('    color: white;')+'\n'
			f.write(line)
			line=str('    text-align: center;')+'\n'
			f.write(line)
			line=str('    padding: 14px 16px;')+'\n'
			f.write(line)
			line=str('    text-decoration: none;')+'\n'
			f.write(line)
			line=str('}')+'\n'

			line=str('li a:hover {')+'\n'
			f.write(line)
			line=str('    background-color:#000;')+'\n'
			f.write(line)
			line=str('}')+'\n'
			f.write(line)
			line=str('</style>')+'\n'
			f.write(line)

			line=str('<ul>')+'\n'
			f.write(line)

				
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/JARVIS.html">HOME</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/JVASP.html">JARVIS-DFT</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/periodic.html">JARVIF-FF</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/documentation.html">DOCUMENTATION</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/others.html">OTHER RESOURCES</a></li>')+'\n'
			f.write(line)
			line=str('	<li ><a href="http://www.ctcms.nist.gov/~knc6/contact.html">CONTACT</a></li>')+'\n'
			f.write(line)
			line=str(' </ul>')+'\n'
			f.write(line)





			line=str("<head>")+'\n'
			f.write(line)
                        finder = SpacegroupAnalyzer(p_strt)
                        ##sgp=str(finder.get_spacegroup_number())
                        sgp=str(finder.get_spacegroup_symbol())
                        sgp.replace('/', 'slash')
                        print ("Checking ITssssssssssssssssssss 1L",file)
                        if '1L' in file:
                               typ=str('1L')
                        elif '2L' in file:
                               typ=str('2L')
                        elif 'Mol' in file:
                               typ=str('Mol')
                        else:
                               typ=str('Bulk')
                        line=str('<h2 style="color:blue;">')+str("Structural formula: ")+str(formula)+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Functional: ")+str(tag)+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Space group number: ")+sgp+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+str("Calculation type: ")+typ+str("</h2>")
			f.write(line)
                        line=str('<h2 style="color:blue;">')+ str(" JARVIS ID: ")+str(name)+str("</h2>")
			f.write(line)
                        if str(final_enp)!='na':
                           line=str('<h2 style="color:blue;">')+str("Final energy/atom: ")+str(final_enp)+str("</h2>")
		    	   f.write(line)
                        path_to_file=str("VASP-FIGS/")+str(name)+str("/")+str("MAIN-RELAX.zip")
                        #line=str('<a href=')+str(path_to_file)+str(' >Download files</a>')
                        line=str('<a href=')+str(path_to_file)+str(' download=')+str(name)+str('>Download files</a>')
			f.write(line)
			line=str("</script>")+'\n'
			f.write(line)
			line=str('<script type="text/javascript"  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>')+'\n'
			f.write(line)
			line=str('<script type="text/javascript" src="JSmol.min.js"></script>')+'\n'
			f.write(line)
			line=str('<script type="text/javascript">')+'\n'
			f.write(line)
			line=str("Info = {")+'\n'
			f.write(line)
			line=str("	width: 400,")+'\n'
			f.write(line)
			line=str("	height: 400,")+'\n'
			f.write(line)
			line=str("	debug: false,")+'\n'
			f.write(line)
			line=str('	color: "0xC0C0C0",')+'\n'
			f.write(line)
			line=str("  disableJ2SLoadMonitor: false,")+'\n'
			f.write(line)
			line=str("  disableInitialConsole: true,")+'\n'
			f.write(line)
			line=str("	addSelectionOptions: false,")+'\n'
			f.write(line)
			line=str('	use: "HTML5",')+'\n'
			f.write(line)
			line=str("	readyFunction: null,")+'\n'
			f.write(line)
			   #line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
			line=str('	script: "load VASP-FIGS/')+str(name)+str("/")+str(name)+str('_prem.cif"')+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("loadJmol = function() {")+'\n'
			f.write(line)
			line=str('  $("#mathdiv").show();')+'\n'
			f.write(line)
			line=str('	$("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))')+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("$(document).ready(function(){")+'\n'
			f.write(line)
			line=str("  checkMathJax(loadJmol);")+'\n'
			f.write(line)
			line=str("});")+'\n'
			f.write(line)
			line=str("checkMathJax = function(fReady) {")+'\n'
			f.write(line)
			line=str('  return (!$("#MathJax_Message")[0] ')+'\n'
			f.write(line)
			line=str("    || $(")+str("'")+str('#MathJax_Message:contains("Loading")')+str("'")+str(")[0]")+'\n'
			f.write(line)
			line=str("    ? setTimeout(function(){checkMathJax(fReady)},10)")+'\n'
			f.write(line)
			line=str("    : fReady());")+'\n'
			f.write(line)
			line=str("}")+'\n'
			f.write(line)
			line=str("</script>")+'\n'
			f.write(line)
			line=str("</head>")+'\n'
			f.write(line)
			line=str("<body>")+'\n'
			f.write(line)
			line=str('<table><tr><td  style="width:450px">')+'\n'
			#line=str('<table><tr><td valign=top style="width:450px">')+'\n'
			#f.write(line)
			#line=str("</td><td>")+'\n'
			f.write(line)
			line=str("<span id=mydiv></span>")+'\n'
			f.write(line)
			line=str("</td></tr></table>")+'\n'
			f.write(line)



			#en_convg=EN_fig(pref='ENCUT',storedir=dir_name)
			#search,k_convg=K_fig(pref='KPOINTS',storedir=dir_name)
			try:
			   group=group+1
			   search,k_convg=K_fig(pref='KPOINTS',storedir=dir_name)

			   en_convg=EN_fig(pref='ENCUT',storedir=dir_name)
			   print ("search===",search,k_convg)
			   line=str(key)+str(',')+str(search)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           info.append(str(search))
                           info.append(str(name))
                           info.append(str(formula))
                           info.append(str(tag))
                           xc_info['search']=str(search)
                           xc_info['name']=str(name)
                           xc_info['formula']=str(formula)
                           xc_info['tag']=str(tag)
			   print ("line======",line)
			   group=group+1
			   line=str(key)+str(',')+str(name)+str(',')+str(period)+str(',')+str(group)+'\n'
			   #df.write(line)
                           ref='NA'
			   ref_f=open("REFERENCE","r")
			   #ref_f=open("REFERENCE","r")
			   for rline in ref_f:
			       r=(str(rline))#.split('\n')[0]
                               if 'mp' in r and 'L' not in r:
                                  try:
                                     ref=(str(r)).split('\n')[0]
                                  except:
                                      pass
                               if '1L' in file: #str(l.split('.vasp')[0]).split('POSCAR-')[1]
                                   #r=str(rline.split('.vasp')[0]).split('POSCAR-')[1]
                                   #import sys
                                   #sys.exit()
                                   try:
                                      r=str(rline.split('.vasp')[0]).split('POSCAR-')[1]
			              #ref=str(str((str(line)).split('\n')[0]).split('.vasp')[0])split('POSCAR-')[1]
                                   except:
                                      pass
                               all_refs.append(str(r))
			   #for line in ref_f:
			   #    ref=(str(line))#.split('\n')[0]
			   #    #ref=(str(line)).split('\n')[0]
                           #    if '1L' in file: #str(l.split('.vasp')[0]).split('POSCAR-')[1]
                           #        try:
                           #           ref=str(line.split('.vasp')[0]).split('POSCAR-')[1]
			   #           #ref=str(str((str(line)).split('\n')[0]).split('.vasp')[0])split('POSCAR-')[1]
                           #        except:
                           #           pass
                           #    all_refs.append[str(ref)]
			   group=group+1
			   line=str(key)+str(',')+str(ref)+str(',')+str(period)+str(',')+str(group)+'\n'
                           #info.append(str(ref))
                           icsd='na'
                           #icsd=get_struct_from_mp(ref)
                           try:
                              icsd=get_struct_from_mp(ref)
                           except:  
                              pass
                           print "icsd2=",icsd
                           if '1L' in file:
                               typ=str('1L')
                           elif '2L' in file:
                               typ=str('2L')
                           elif 'Mol' in file:
                               typ=str('Mol')
                           else:
                               typ=str('Bulk')
                               if 'mp' in ref:
                                 try:
                                    icsd=get_struct_from_mp(ref)
                                    print 'icsd=',icsd
                                 except:
                                    pass
                           info.append(typ)
                           xc_info['typ']=str(typ)
                           xc_info['bandgap']=str('na')
                           xc_info['bandgapHSE']=str('na')
                           xc_info['BV']=str('na')
                           xc_info['GV']=str('na')
                           xc_info['SG']=str(sgp)#str('na')
                           xc_info['CWD']=str(os.getcwd())
                           xc_info['enp']=str('na')
                           xc_info['dif']=str('na')
                           xc_info['np']=str('na')
                           N_A=str('na')
                           info.append(N_A)
                           info.append(N_A)
                           info.append(N_A)
                           info.append(N_A)
                           try:
                              info.append(sgp)
                           except:
                              info.append(N_A)
                              pass
                        
                           info.append(N_A)
                           #info.append(N_A)
#["Search key","JVASP ID","Formula","Functional","REFID","Type","Bandgap&nature","Bv (GPa)","Gv (GPa)","SG","EN"]
			   #df.write(line)

			   line=str('<!--Source:')+str(os.getcwd())+str(' -->')+'\n'
			   f.write(line)
			   line=str('<h2 style="color:blue;">Convergence</h2>')+'\n'
			   f.write(line)
                           #line=str('<p > This is convergenece test </p>')+'\n'
                           line=str('<p > Calculations are done using VASP software. Convergence on KPOINTS and ENCUT is done with respect to total energy of the system within 0.001 eV tolerance. Please note convergence on KPOINTS and ENCUT is generally done for target properties, but here we assume energy-convergence with 0.001 eV should be sufficient for other properties also. The points on the curves are obtained with single-point calculation (nuber of ionic steps, NSW=1). However, for very accurate calculations, NSW>1 might be needed.   </p>')+'\n'
	                   f.write(line)
			   line=str('    <div class="row">')+'\n'
			   f.write(line)
			   line=str('      <div class="col s12 m6">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/KDen.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('        <image src="VASP-FIGS/')+str(name)+str('/Encut.png" width="450" height="350">')+'\n'
			   f.write(line)
			   line=str('      </div>')+'\n'
			   f.write(line)
			except:
                              pass
        os.chdir(cwd)
        group=0
        period=0
        for el in info:    
           group=group+1
           line=str(key)+str(',')+str(el)+str(',')+str(period)+str(',')+str(group)+'\n'
           #df.write(line)
           		    
        if info !=[]:

           Infarray.append(info)


   print  "infarray",Infarray
   #sort_infarry=[]
   #sort_infarry=sorted(Infarray,key=lambda student: student[0])
   #sort_infarry=(Infarray)



##   fileXC=str('XCinfo32d')+'_data.json'
##   fjXC=open(fileXC,'w')

##   fjXC.write(json.dumps(xc_info_arr,cls=MontyEncoder,indent=4))
##   fjXC.close()


   sort_infarry=sorted(Infarray)
   print  "sorted_infarray", sort_infarry
   keys=[]
   new=0
   count1=0

   for i in range(0,len(sort_infarry)):
       if sort_infarry[i][0] not in keys:
           keys.append(sort_infarry[i][0])
           new=new+1
           count1=0
       elif sort_infarry[i][0]  in keys:
           count1=count1+1
       for j in range(0,len(sort_infarry[i])):
           #print Intfarray[i][j],count+1,j+1
           line= str(new)+str(',')+str(sort_infarry[i][j])+str(',')+str(count1+1)+str(',')+str(j+1)+'\n'
           df.write(line)


#   for i,el1 in enumerate(Infarray):

#   for i in range(0,len(sort_infarry)):
#       if el[0] not in keys:
#           keys.append(el[0])
#           new=new+1
#           count1=0
#       elif el[0]  in keys:
#           count1=count1+1
#       for j,el2 in enumerate(el1):
#           #print Intfarray[i][j],count+1,j+1
#           line= str(new)+str(',')+str(el2)+str(',')+str(count1+1)+str(',')+str(j+1)+'\n'
#           df.write(line)










   df.close()
       #main(chdir=file,storedir=dir_name)
#Ramanplot('/data/knc6/NIST2/For2D/MoTe2_mp-602/RAMANDIR-bulk@mp_602/vasp_raman.dat', '/data/knc6/NIST2/For2D/STORE/JV2D-1')
#icsd=get_struct_from_mp('mp-134')
#print icsd
loop_over() 


periodic(html='jjaa.html',csv='DFT_dathold32d.csv',header=["Search key","JVASP ID","Formula","Functional","Type","Bandgap&nature","Bv (GPa)","Gv (GPa)","SG","EN"])

#>>> for el in ss:
#...    if len(str(el[0]).split('-'))>1 and el[0] not in key2:
#...        key2.append(el[0])
#...    elif '-' not in el[0] and el[0] not in key1:
#...        key1.append(el[0])
#...
#>>> print key1

#periodic(html='JARVIS-2D.html',csv='DFT.csv')
