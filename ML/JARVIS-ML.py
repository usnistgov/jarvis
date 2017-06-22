# filter warnings messages from the notebook
import warnings
warnings.filterwarnings('ignore')
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import pickle,os

from matminer.descriptors.composition_features import get_pymatgen_descriptor
from pymatgen import Composition

descriptors = ['row', 'group', 'atomic_mass', 'atomic_radius', 'boiling_point', 'melting_point', 'X']
import glob

#read_files = glob.glob("*.json")
#with open("merged_file.json", "wb") as outfile:
#    outfile.write('[{}]'.format(
#        ','.join([open(f, "rb").read() for f in read_files])))

#data = loadfn('all_info_data.json', cls=MontyDecoder)

def prepare_x(struct=None):
    val=[]
    for i in descriptors:
        d=np.mean(get_pymatgen_descriptor(struct.composition,i))
        val.append(d)
    return val

def train_data():
        data = loadfn('all_info_0_data.json', cls=MontyDecoder)
	X=[]
	Y=[]
	for i in data:
	    #print i['final_str']
	    #print i['opt_info']['refrx'][0]
	    #print i['opt_info']['refry'][0]
	    #print i['opt_info']['refrz'][0]
	    #print i['bbandgap']
	    #print i['elastic']
            try:
	       val=prepare_x(i['final_str'])
               temp=[i['opt_info']['refrx'][0],float(i['bbandgap'].split()[0]),i['elastic'][0][0]]
               typ=str(np.array(list(temp)).dtype)
               print typ
               #add=1
               #for j in temp:
               #    if str(j)=='na':
               #        print 'na'
               #        add=0
               #        #print j
               #if add==1:
               if typ=='float64':
               
	               X.append(val)
	               Y.append([i['opt_info']['refrx'][0],float(i['bbandgap'].split()[0]),i['elastic'][0][0]])
            except:
                    pass
	    #print np.mean(get_pymatgen_descriptor(i['final_str'].composition,'X'))
	    #print val
	    #print np.mean(get_pymatgen_descriptor(i['final_str'].composition,'X'))
	    #print str(i['final_str'].composition).map(lambda x: np.mean(get_pymatgen_descriptor(Composition(x), d)))
	#print X,len(X)
	#print Y,len(Y)
	linear_regression = LinearRegression()
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
	linear_regression.fit(X, Y)
        pickle.dump( linear_regression, open( "save.p", "wb" ) )

file=str(os.getcwd())+str('/save.p')
if not os.path.isfile(file):
   train_data()
#train_data()
object_file =pickle.load( open( "save.p", "rb" ) )
#entry='FeOCH'
entry = raw_input('Enter a formula: ')
print "You entered",entry
comp_ent=Composition(entry)
ent_X=[]
for i in descriptors:
    d=np.mean(get_pymatgen_descriptor(comp_ent,i))
    ent_X.append(d)
#print comp_ent,ent_X
pred=object_file.predict(ent_X)
print 'predicted refr_x',pred[0][0]
print 'predicted bandgap',pred[0][1]
print 'predicted C11',pred[0][2]

