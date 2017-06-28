from pymatgen import Composition, Element
from numpy import zeros, mean

# Training file containing band gaps extracted from Materials Project
# created in previous blog post and linked here
trainFile = open("bandgapDFT.csv","r").readlines()

# Input: pymatgen Composition object
# Output: length-100 vector representing any chemical formula

def naiveVectorize(composition):
       vector = zeros((MAX_Z))
       for element in composition:
               fraction = composition.get_atomic_fraction(element)
               vector[element.Z - 1] = fraction
       return(vector)

# Extract materials and band gaps into lists, and construct naive feature set
materials = []
bandgaps = []
naiveFeatures = []

MAX_Z = 100 # maximum length of vector to hold naive feature set

for line in trainFile:
       split = str.split(line, ',')
       material = Composition(split[0])
       materials.append(material) #store chemical formulas
       naiveFeatures.append(naiveVectorize(material)) #create features from chemical formula
       bandgaps.append(float(split[1])) #store numerical values of band gaps
print materials[0]
print
print bandgaps[0]
print
print naiveFeatures[0]
