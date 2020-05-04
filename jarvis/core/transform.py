from monty.serialization import loadfn
d=loadfn('periodic_table.json')
dd=loadfn('Elements.json')
allkeys=dd['Pr'].keys()
gg=dd

no_there=[]
for i in d.keys():
 if i not in dd.keys():
    no_there.append(i)


for i in no_there:
   gg.setdefault(i,{})

for i in no_there:
  for j in allkeys:
    gg[i][j]=-9999
    for k in d[i].keys():
       if d[i][k]=='no data':
         d[i][k]=-9999
    gg[i]['Z']=d[i]['Atomic no']
    gg[i]['atom_mas']=d[i]['Atomic mass']
    gg[i]['atom_rad']=d[i]['Atomic radius calculated']
    gg[i]['X']=d[i]['X']


import json
f=open('Elements_new.json','w')
f.write(json.dumps(gg))
f.close()




