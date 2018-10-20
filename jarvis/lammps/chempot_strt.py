from __future__ import  unicode_literals, print_function
from pymatgen.ext.matproj import MPRester
import json,os,operator
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.core.periodic_table import get_el_sp, Element

def get_struct_from_mp(formula, MAPI_KEY="", all_structs=False):

    if not MAPI_KEY:
        MAPI_KEY = os.environ.get("MAPI_KEY", "")
        if not MAPI_KEY:
            print('API key not provided')
            print('get API KEY from materialsproject and set it to the MAPI_KEY environment variable. aborting ... ')
            sys.exit()
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        structures = []
        x = {}
        print("\nnumber of structures matching the chemical formula {0} = {1}".format(formula, len(data)))
        print("The one with the the lowest energy above the hull is returned, unless all_structs is set to True")
        for d in data:
            mpid = str(d['material_id'])
            x[mpid] = d['e_above_hull']
            if all_structs:
                structure = m.get_structure_by_material_id(mpid)
                structure.sort()
                structures.append(structure)
        if all_structs:
            return structures
        else:
            key = sorted(x.items(), key=operator.itemgetter(1))[0][0]
            print("The id of the material corresponding to the lowest energy above the hull = {0}".format(key))
            if key:
                return key,m.get_structure_by_material_id(key)
            else:
                return None
if __name__=='__main__':
 f=open('chempot_strt.json','w')
 mem=[]
 for i in Element:
  try:
   key,strt=get_struct_from_mp(i.symbol)
   info={}
   info['element']=str(i.symbol)
   info['mpid']=key
   info['structure']=strt
   mem.append(info)
  except:
   pass
 f.write(json.dumps(mem,cls=MontyEncoder))
 f.close()
     
