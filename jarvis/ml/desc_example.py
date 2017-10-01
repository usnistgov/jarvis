from pymatgen.core.structure import Structure
from get_desc import get_comp_descp
s=Structure.from_file('POSCAR')
X=get_comp_descp(s)
print 'X',X,len(X)

