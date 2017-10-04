import itertools
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from get_des import get_comp_descp
import sys
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from sklearn.feature_selection import SelectKBest
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
import numpy as np
import pickle
def pd_plot(system='Ni-Nb'):


  object_file =pickle.load( open( 'form_en_gbsave_v25h', "rb" ))
  output=[]
  l=system.split('-')
  comb = []
  for i in range(len(l)):
      comb += itertools.combinations(l,i+1)
  comb_list = [ list(t) for t in comb ]
#  comb_list=['Zn','O','Zn-O']
  with MPRester("") as m:
    for i in    comb_list:
        dd='-'.join(i)
        print dd
        data = m.get_data(dd)
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            X=get_comp_descp(struct=structure)
            pred=object_file.predict(X)
            print structure.composition.reduced_formula,pred[0],d['formation_energy_per_atom'],str(d['material_id'])
            output.append(PDEntry(Composition(structure.composition),float(pred[0])))
  pd = PhaseDiagram(output)
  print output
  plotter = PDPlotter(pd, show_unstable=True)
  name=str(system)+str('_phasediagram.png')
  plotter.write_image(name,image_format="png")


#pd_plot()
