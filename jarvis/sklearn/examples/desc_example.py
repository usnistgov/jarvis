# coding: utf-8
from jarvis.sklearn.get_desc import get_comp_descp
from pymatgen.core.structure import Structure
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import numpy as np
import sys


#STEP-1: Check if the array length for POSCAR is 1866
s=Structure.from_file('POSCAR')
X=get_comp_descp(s)
print (X,len(X))

#STEP-2: Download JARVIS-DFT database for bulk materials
# wget --no-check-certificate https://www.ctcms.nist.gov/~knc6/jdft_3d.json.tgz
#tar -xvzf jdft_3d.json.tgz
#jdft_3d-11-11-2017.json

#STEP-3: load data 
data=loadfn('jdft_3d-11-11-2017.json',cls=MontyDecoder)
print ('data length',len(data)) #12888
print ('available keys',data[0].keys())

def isfloat(value):
  #Simple function to check if the data is available/float
  try:
    float(value)
    return True
  except ValueError:

#STEP-4:  Let's train a Formation energy ML

    return False
#for i in data:
#  print (i['form_enp'],i['mpid'],i['jid'],i['final_str'],type(i['final_str']))
#  sys.exit()

#Making descriptors and target data exportable to sklearn
X=[]
Y=[]
for ii,i in enumerate(data):
 if ii<5:
  y=i['form_enp']
  if isfloat(y):
    x=get_comp_descp(i['final_str'])
    X.append(x)
    Y.append(y)

#STEP-5: Fit and predict
X=np.array(X).astype(np.float)
Y=np.array(Y).astype(np.float)
est= GradientBoostingRegressor()
pipe=Pipeline([ ("fs", VarianceThreshold()),("est", est)])
pipe.fit(X,Y)

test_x=get_comp_descp(data[-1]['final_str'])
print (pipe)
print (pipe.predict([test_x]),data[-1]['form_enp'])
