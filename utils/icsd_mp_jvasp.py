from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.json import MontyEncoder, MontyDecoder
import json
f=open('JAVSP-MP-latt.json','r')
data=json.load(f,cls=MontyDecoder)
f.close()
from monty.serialization import loadfn, dumpfn
datamp = loadfn('/home/knc6/bin/MPall_datacopy.json', cls=MontyDecoder)

#https://github.com/hswayne77/python-mae/blob/master/mae.py
def mae(pred_list, true_list):
    if len(pred_list) != len(true_list):
        raise Exception('Error: number of elements not match!')
    return sum(map(lambda t:float(t[0]-t[1]),zip(pred_list, true_list)))/len(true_list)
#for i in data:
#   print i['JVASP'][0] 
icsd_a_x=[]
icsd_a_y=[]
icsd_b_x=[]
icsd_b_y=[]
icsd_c_x=[]
icsd_c_y=[]
j_a_x=[]
j_a_y=[]
j_b_x=[]
j_b_y=[]
j_c_x=[]
j_c_y=[]
m_a_x=[]
m_a_y=[]
m_b_x=[]
m_b_y=[]
m_c_x=[]
m_c_y=[]


m_a_err=[]
j_a_err=[]
m_b_err=[]
j_b_err=[]
m_c_err=[]
j_c_err=[]
m_a_id=[]
m_b_id=[]
m_c_id=[]
avoid=['mp-8677','mp-21405','mp-19846','mp-27734','mp-6462','mp-19755','mp-23007','mp-6462','mp-24814','mp-383','mp-1797','mp-15679','mp-19079','mp-23294','mp-23312']
for d in datamp:
  mmpid= str(d['mp_id'])
  icsd=d["icsd"]
  if icsd!=None and icsd!=[] :
    for i in data:
       jmpid= i["mpid"]
       if mmpid==jmpid and mmpid not in avoid:
         ini= (d['ini_structure']).get_primitive_structure()
         ini_cvn=SpacegroupAnalyzer(ini).get_conventional_standard_structure()
         if len(ini_cvn)==len(ini):
             ini=ini_cvn
         ia,ib,ic=ini.lattice.abc
         ja=i['JVASP'][0]
         jb=i['JVASP'][1]
         jc=i['JVASP'][2]
         mpa=i['MP'][0]
         mpb=i['MP'][1]
         mpc=i['MP'][2]
         m_a=abs(100*float(ia-mpa)/float(ia))
         if m_a>=5:
			 icsd_a_x.append(ia)
			 icsd_a_y.append(ia)
			 j_a_x.append(ia)
			 j_a_y.append(ja)
			 m_a_x.append(ia)
			 m_a_y.append(mpa)
			 m_a_id.append(mmpid)
                         if m_a>35:
                            print mmpid,'error in a gt 35',m_a,ia,mpa
		 
			 j_a=abs(100*float(ia-ja)/float(ia))
			 m_a_err.append(m_a)
			 j_a_err.append(j_a)
                         if mmpid=='mp-8677':
                              print "icsd",ia,ib,ic,"MP",mpa,mpb,mpc,"J",ja,jb,jc
         m_b=abs(100*float(ib-mpb)/float(ib))
         if m_b>=5:

			 icsd_b_x.append(ib)
			 icsd_b_y.append(ib)
			 j_b_x.append(ib)
			 j_b_y.append(jb)
			 m_b_x.append(ib)
			 m_b_y.append(mpb)

                         if m_b>35:
                            print mmpid,'error in b gt 35',m_b,ib,mpb
			 j_b=abs(100*float(ib-jb)/float(ib))
			 m_b_err.append(m_b)
			 j_b_err.append(j_b)
			 m_b_id.append(mmpid)


         m_c=abs(100*float(ic-mpc)/float(ic))
         if m_c>=5:
			 icsd_c_x.append(ic)
			 icsd_c_y.append(ic)
			 j_c_x.append(ic)
			 j_c_y.append(jc)
			 m_c_x.append(ic)
			 m_c_y.append(mpc)
                         if m_c>35:
                            print mmpid,'error in c gt 35',m_c,ic,mpc

			 j_c=abs(100*float(ic-jc)/float(ic))
			 m_c_err.append(m_c)
			 j_c_err.append(j_c)
			 m_c_id.append(mmpid)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
font = {'size'   : 18}
import matplotlib
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from pymatgen.util.plotting_utils import get_publication_quality_plot
from matplotlib.gridspec import GridSpec

the_grid = GridSpec(2, 3)
plt = get_publication_quality_plot(18, 14)
plt.subplot(the_grid[0, 0])
plt.plot(icsd_a_x,icsd_a_y,'b-')
plt.plot(j_a_x,j_a_y,'*')
plt.plot(m_a_x,m_a_y,'.')

plt.subplot(the_grid[0, 1])
plt.plot(icsd_b_x,icsd_b_y,'b-')
plt.plot(j_b_x,j_b_y,'*')
plt.plot(m_b_x,m_b_y,'.')

plt.subplot(the_grid[0, 2])
#plt.subplot(the_grid[0, 2], aspect=1)
plt.plot(icsd_c_x,icsd_c_y,'b-')
plt.plot(j_c_x,j_c_y,'*')
plt.plot(m_c_x,m_c_y,'.')
plt.savefig('C_COMP.png')
plt.close()


the_grid = GridSpec(4, 3)
the_grid.update( hspace=0.75)
the_grid.update( wspace=0.4)
plt = get_publication_quality_plot(18, 14)

plt.subplot(the_grid[0, 0])
plt.plot(icsd_a_x,icsd_a_y,'b-')
plt.plot(m_a_x,m_a_y,'.')
#for i,j,k,l in zip(m_a_err,m_a_x,m_a_y,m_a_id):
#     if i>25:
#         plt.text(j,k,l)
mae_m_a=round(mae(m_a_y,m_a_x),2)
ttl=str('(a)')+str(' \n MAE:')+str(mae_m_a)
plt.title(ttl)
plt.xlabel(r'a (ICSD)($\AA$)')
plt.ylabel(r'a (MP)($\AA$)')

plt.subplot(the_grid[1, 0])
plt.bar(icsd_a_x,m_a_err,width=.05,color='r')
plt.ylim([0,100])
plt.title('(b)')
plt.xlabel(r'a (ICSD)($\AA$)')
plt.ylabel('MP-Rel. perc. error')

plt.subplot(the_grid[2, 0])
plt.plot(icsd_a_x,icsd_a_y,'b-')
plt.plot(j_a_x,j_a_y,'*')
mae_j_a=round(mae(j_a_y,j_a_x),2)
ttl=str('(c)')+str(' \n MAE:')+str(mae_j_a)
plt.title(ttl)
plt.xlabel(r'a (ICSD)($\AA$)')
plt.ylabel(r'a (JVASP)($\AA$)')


plt.subplot(the_grid[3, 0])
plt.bar(icsd_a_x,j_a_err,width=.05,color='b')
plt.ylim([0,100])
plt.xlabel(r'a (ICSD)($\AA$)')
plt.ylabel('JVASP-Rel. perc. error')
plt.title('(d)')

#####################################

plt.subplot(the_grid[0, 1])
plt.plot(icsd_b_x,icsd_b_y,'b-')
plt.plot(m_b_x,m_b_y,'.')
#for i,j,k,l in zip(m_b_err,m_b_x,m_b_y,m_b_id):
#     if i>25:
#         plt.text(j,k,l)
mae_m_b=round(mae(m_b_y,m_b_x),2)
ttl=str('(e)')+str(' \n MAE:')+str(mae_m_b)
plt.title(ttl)

plt.xlabel(r'b (ICSD)($\AA$)')
plt.ylabel(r'b (MP)($\AA$)')


plt.subplot(the_grid[1, 1])
plt.bar(icsd_b_x,m_b_err,width=.05,color='r')
plt.ylim([0,100])
plt.title('(f)')
plt.xlabel(r'b (ICSD)($\AA$)')
plt.ylabel('MP-Rel. perc. error')

plt.subplot(the_grid[2, 1])
plt.plot(icsd_b_x,icsd_b_y,'b-')
plt.plot(j_b_x,j_b_y,'*')
mae_j_b=round(mae(j_b_y,j_b_x),2)
ttl=str('(g)')+str(' \n MAE:')+str(mae_j_b)
plt.title(ttl)
plt.xlabel(r'b (ICSD)($\AA$)')
plt.ylabel(r'b (JVASP)($\AA$)')


plt.subplot(the_grid[3, 1])
plt.bar(icsd_b_x,j_b_err,width=.05,color='b')
plt.ylim([0,100])
plt.xlabel(r'b (ICSD)($\AA$)')
plt.ylabel('JVASP-Rel. perc. error')
plt.title('(d)')
plt.title('(h)')

#####################################


plt.subplot(the_grid[0, 2])
plt.plot(icsd_c_x,icsd_c_y,'b-')
plt.plot(m_c_x,m_c_y,'.')
#for i,j,k,l in zip(m_c_err,m_c_x,m_c_y,m_c_id):
#     if i>25:
#         plt.text(j,k,l)
mae_m_c=round(mae(m_c_y,m_c_x),2)
ttl=str('(i)')+str(' \n MAE:')+str(mae_m_c)
plt.title(ttl)
plt.xlabel(r'c (ICSD)($\AA$)')
plt.ylabel(r'c (MP)($\AA$)')

plt.subplot(the_grid[1, 2])
plt.bar(icsd_c_x,m_c_err,width=.05,color='r')
plt.ylim([0,100])
plt.title('(j)')
plt.xlabel(r'c (ICSD)($\AA$)')
plt.ylabel('MP-Rel. perc. error')

plt.subplot(the_grid[2, 2])
plt.plot(icsd_c_x,icsd_c_y,'b-')
plt.plot(j_c_x,j_c_y,'*')
mae_j_c=round(mae(j_c_y,j_c_x),2)
ttl=str('(k)')+str(' \n MAE:')+str(mae_j_c)
plt.title(ttl)
plt.xlabel(r'c (ICSD)($\AA$)')
plt.ylabel(r'c (JVASP)($\AA$)')

plt.subplot(the_grid[3, 2])
plt.bar(icsd_c_x,j_c_err,width=.05,color='b')
plt.title('(l)')
plt.xlabel(r'c (ICSD)($\AA$)')
plt.ylabel('JVASP-Rel. perc. error')
plt.ylim([0,100])


#####################################
plt.savefig('C_COMP2.png')
plt.close()
###############
############

the_grid = GridSpec(2, 3)
plt = get_publication_quality_plot(18, 14)

plt.subplot(the_grid[0, 0])
plt.plot(icsd_a_x,icsd_a_y,'b-')
plt.plot(m_a_x,m_a_y,'.')
plt.xlabel(r'a(ICSD)($\AA$)')
plt.ylabel(r'a(MP)($\AA$)')
plt.title('(a)')

plt.subplot(the_grid[1, 0])
plt.plot(icsd_a_x,icsd_a_y,'b-')
plt.plot(j_a_x,j_a_y,'.')
plt.xlabel(r'a(ICSD)($\AA$)')
plt.ylabel(r'a(JVASP)($\AA$)')
plt.title('(d)')
for i,j,k,l in zip(j_a_err,j_a_x,j_a_y,m_a_id):
     if i>20 and j>15:
         plt.text(j,k,l)



#####################################

plt.subplot(the_grid[0, 1])
plt.plot(icsd_b_x,icsd_b_y,'b-')
plt.plot(m_b_x,m_b_y,'.')
plt.xlabel(r'b(ICSD)($\AA$)')
plt.ylabel(r'b(MP)($\AA$)')
plt.title('(b)')


plt.subplot(the_grid[1, 1])
plt.plot(icsd_b_x,icsd_b_y,'b-')
plt.plot(j_b_x,j_b_y,'.')
plt.xlabel(r'b(ICSD)($\AA$)')
plt.ylabel(r'b(JVASP)($\AA$)')
plt.title('(e)')

for i,j,k,l in zip(j_b_err,j_b_x,j_b_y,m_b_id):
     if i>20 and j>15:
         plt.text(j,k,l)


#####################################


plt.subplot(the_grid[0, 2])
plt.plot(icsd_c_x,icsd_c_y,'b-')
plt.plot(m_c_x,m_c_y,'.')
plt.xlabel(r'c(ICSD)($\AA$)')
plt.ylabel(r'c(MP)($\AA$)')
plt.title('(c)')


plt.subplot(the_grid[1, 2])
plt.plot(icsd_c_x,icsd_c_y,'b-')
plt.plot(j_c_x,j_c_y,'.')
plt.xlabel(r'c(ICSD)($\AA$)')
plt.ylabel(r'c(JVASP)($\AA$)')
plt.title('(f)')

for i,j,k,l in zip(j_c_err,j_c_x,j_c_y,m_c_id):
     if i>20 and j>15:
         plt.text(j,k,l)


#####################################
plt.savefig('C_COMP3.png')
plt.close()

print "MP A error",max(m_a_err),min(m_a_err),max(j_a_err),min(j_a_err)
print "MP B error",max(m_b_err),min(m_b_err),max(j_b_err),min(j_b_err)
print "MP C error",max(m_c_err),min(m_c_err),max(j_c_err),min(j_c_err)
print "MP A error",max(m_a_err),min(m_a_err),max(j_a_err),min(j_a_err)
print "MP B error",max(m_b_err),min(m_b_err),max(j_b_err),min(j_b_err)
print "MP C error",max(m_c_err),min(m_c_err),max(j_c_err),min(j_c_err)
