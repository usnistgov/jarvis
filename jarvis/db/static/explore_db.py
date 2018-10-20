from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn

"""
 Download 3D material data from https://figshare.com/articles/jdft_3d-7-7-2018_json/6815699
 Download 2D material data from https://figshare.com/articles/jdft_2d-7-7-2018_json/6815705 

 Website: https://jarvis.nist.gov
 https://www.nature.com/articles/s41598-017-05402-0
 https://www.nature.com/articles/sdata201882
 https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.014107
"""


if __name__=='__main__':
 d=loadfn('jdft_3d.json',cls=MontyDecoder)
 #d=loadfn('jdft_2d.json',cls=MontyDecoder)

 print ("Total size",len(d))
 for i in d:
    print ("JARVIS id", i['jid'],"link",str("https://www.ctcms.nist.gov/~knc6/jsmol/")+str(i['jid'])+str(".html"))
    print ("Materials project id"), i['mpid'],"link",str("https://www.materialsproject.org/materials/")+str(i['mpid'])
    print ("final energy eV",i['fin_en'])
    print ("magnetic moment, Bohr Mag.",i['magmom'])
    print ("optb88vdw bandgap eV", i["op_gap"])
    print ("Tran-Blaha MBJ bandgap  eV",i["mbj_gap"])
    print ("static dielectric function in x,optb88vdw", i["epsx"]) #Careful about 2D
    print ("static dielectric function in x,MBJ", i["mepsx"]) #Careful about 2D
    print ("bulk modulus GPa", i["kv"]) #Careful about 2D
    print ("Shear modulus GPa", i["gv"]) #Careful about 2D
    break
