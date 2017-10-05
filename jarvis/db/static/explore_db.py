from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn

#https://www.ctcms.nist.gov/~knc6/jdft_3d.json.tgz
#https://www.ctcms.nist.gov/~knc6/jdft_2d.json.tgz
#tar -xvzf *.tgz


d=loadfn('jdft_3d.json',cls=MontyDecoder)
#d=loadfn('jdft_2d.json',cls=MontyDecoder)

print "Total size",len(d)
for i in d:
    print "JARVIS id", i['jid'],"link",str("https://www.ctcms.nist.gov/~knc6/jsmol/")+str(i['jid'])+str(".html")
    print "Materials project id", i['mpid'],"link",str("https://www.materialsproject.org/materials/")+str(i['mpid'])
    print  "final energy eV",i['fin_en']
    print "optb88vdw bandgap eV", i["op_gap"]
    print "Tran-Blaha MBJ bandgap  eV",i["mbj_gap"]
    print "static dielectric function in x,optb88vdw", i["epsx"] #Careful about 2D
    print "static dielectric function in x,MBJ", i["mepsx"] #Careful about 2D
    print "bulk modulus GPa", i["kv"] #Careful about 2D
    print "Shear modulus GPa", i["gv"] #Careful about 2D
    break
