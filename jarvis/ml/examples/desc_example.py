# coding: utf-8
from pymatgen.core.structure import Structure
from jarvis.ml.get_desc import get_comp_descp
import pickle
s=Structure.from_file('POSCAR')
X=get_comp_descp(s)
with open('form_en_gbsave_v25h', 'rb') as f:
    data = pickle.load(f)
f=open('form_en_gbsave_v25h', 'rb')
#d=pickle.load(f, encoding='utf-8')
f.close()
with open('form_en_gbsave_v25h','rb') as f:
   object_file=pickle.load(f)
#object_file =pickle.load( open( 'form_en_gbsave_v25h', "rb" ))
print ('X',X,len(X))
p=object_file.predict(X)[0]
print ('formation energy/atom prediction',p,' eV')

