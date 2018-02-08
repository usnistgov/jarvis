import math,itertools
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

class Pos(object):
  def __init__(self,frac_coords=None,cart_coords=None,
    elements=None,box=None,pbc=None,
    charge=None,props=None,comment='jarvis'):
    self.frac_coords=frac_coords 
    self.cart_coords=cart_coords 
    self.box=box
    self.charge=charge
    self.elements=elements
    self.pbc=pbc
    self.props=props
    self.comment=comment
    if self.frac_coords is not None and self.box is not None and self.cart_coords is None:
         self.cart_coords=self.get_cart(self.frac_coords)
    if self.cart_coords is not None and self.box is not None and self.frac_coords is None:
         self.frac_coords=self.get_frac(self.frac_coords)
  def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}={value}\n".format(key=key, value=self.__dict__[key]))
 
        return ''.join(sb)
 
  def read_poscar(self,fname):
     f=open(fname,'r')
     lines=f.read().splitlines()
     f.close()
     #print lines
     box=np.zeros((3,3))
     for i,ii in enumerate(lines):
       if i==0:
        comm=ii.split('\n')[0]
       elif i==1:
         scal=float(ii.split()[0])
       elif i==2:
          tmp=ii.split()
          box[0][0]=tmp[0]
          box[0][1]=tmp[1]
          box[0][2]=tmp[2]
       elif i==3:
         tmp=ii.split()
         box[1][0]=tmp[0]
         box[1][1]=tmp[1]
         box[1][2]=tmp[2]
       elif i==4:
         tmp=ii.split()
         box[2][0]=tmp[0]
         box[2][1]=tmp[1]
         box[2][2]=tmp[2]
       elif i==5:
         symbs=str(ii).split()
         box=np.array(box)*scal
       elif i==6:
          types=ii.split()
          types=np.array([ int(i) for i in types])
          nat=np.sum(types)
          pos=np.zeros((nat,3))
       elif i==7:
          cart_frac=ii.split('\n')
       else:
          pos[i-8][0]=float(ii.split()[0]) 
          pos[i-8][1]=float(ii.split()[1]) 
          pos[i-8][2]=float(ii.split()[2]) 
     all_symbs=[]
     for i,ii in enumerate(types):
       for j in range(ii):
         all_symbs.append(symbs[i])
     
     self.box=box
     type=cart_frac[0]
     if type.startswith('d') or type.startswith('D'): 
       frac_coords=pos
       cart_coords=self.get_cart(pos)
     if type.startswith('c') or type.startswith('C'): 
       cart_coords=pos
       frac_coords=self.get_frac(pos)
    
     return Pos(frac_coords=frac_coords,cart_coords=cart_coords,comment=comm,elements=(all_symbs),box=box)
         
       
    
  def supercell(self,dim=[]):
    dim=np.array(dim)
    coords=self.frac_coords
    box=self.box
    all_symbs=self.elements
    #print 'before super',all_symbs
    nat=len(coords)
    new_nat=nat*dim[0]*dim[1]*dim[2]
    new_coords=np.zeros((new_nat,3))
    new_symbs= [] #np.chararray((new_nat))
    
    count=0
    for i in range(nat):
      for j in range(dim[0]):
       for k in range(dim[1]):
        for l in range(dim[2]):
          new_coords[count][0]=coords[i][0]+j    
          new_coords[count][1]=coords[i][1]+k    
          new_coords[count][2]=coords[i][2]+l 
          #new_symbs[count]=all_symbs[i]
          new_symbs.append(all_symbs[i])
          count=count+1
    box=np.multiply(box,dim)  
    cart_coords=self.get_cart(new_coords)
    return Pos(frac_coords=new_coords,elements=new_symbs,box=box)

           
      

     
  def get_cart(self,coords):
      new_coords=np.zeros(coords.shape)
     
      lat=self.box
      for i in range(coords.shape[0]):
       new_coords[i][0]=lat[0][0]*coords[i][0]+lat[1][0]*coords[i][1]+lat[2][0]*coords[i][2]
       new_coords[i][1]=lat[0][1]*coords[i][0]+lat[1][1]*coords[i][1]+lat[2][1]*coords[i][2]
       new_coords[i][2]=lat[0][2]*coords[i][0]+lat[1][2]*coords[i][1]+lat[2][2]*coords[i][2]
      return new_coords
 

  def get_frac(self,coords):
      new_coords=np.zeros(coords.shape)
      coords=self.cart_coords
      lat=self.box
      lat=np.linalg.inv(self.box)
      for i in range(len(coords)):
       new_coords[i][0]=lat[0][0]*coords[i][0]+lat[1][0]*coords[i][1]+lat[2][0]*coords[i][2]
       new_coords[i][1]=lat[0][1]*coords[i][0]+lat[1][1]*coords[i][1]+lat[2][1]*coords[i][2]
       new_coords[i][2]=lat[0][2]*coords[i][0]+lat[1][2]*coords[i][1]+lat[2][2]*coords[i][2]
      return new_coords



  def write_poscar(self,type='cart',fname=''):
    f=open(fname,'w')
    f.write('%s\n' %self.comment)
    f.write('%s\n' %'   1.0')
    box=self.box
    elements=self.elements
    if type.startswith('c') or type.startswith('C'):
       coords=self.cart_coords
    if type.startswith('d') or type.startswith('D'):
       coords=self.frac_coords
    
    f.write('%12.6f%12.6f%12.6f \n' %(box[0][0],box[0][1],box[0][2]))
    f.write('%12.6f%12.6f%12.6f \n' %(box[1][0],box[1][1],box[1][2]))
    f.write('%12.6f%12.6f%12.6f \n' %(box[2][0],box[2][1],box[2][2]))
    
    order = np.argsort(elements)
    els=elements[order]
    cords=coords[order]
    unique_items, counts = np.unique(els, return_counts=True)
    f.write('%5s\n' %(' '.join(unique_items)))
    f.write('%5s\n' %(' '.join(np.array(counts,dtype='string'))))
    if type.startswith('d') or type.startswith('D'):
      f.write('%s\n' %'direct')
    if type.startswith('c') or type.startswith('C'):
      f.write('%s\n' %'cartesian')
    for i,j in zip(coords,els):
      f.write('%12.6f%12.6f%12.6f %s \n' %(i[0],i[1],i[2],j))
     
    f.close()

  def verlet(self,dim=[4,4,4]):
     pp=self.elements
     znm=0
     deg_arr=[]
     bond_arr=[]
     lat=self.box
     
     super=self.supercell(dim=dim)
     coords=super.frac_coords
     all_symbs=super.elements
     dim05=[float(i/2.) for i in dim]
     nat=len(coords)
     #print nat
     nn=np.zeros((nat),dtype='int')
     max_n=450 #maximum number of neighbors
     dist=np.zeros((max_n,nat))
     nn_id=np.zeros((max_n,nat),dtype='int')
     bondx=np.zeros((max_n,nat))
     bondy=np.zeros((max_n,nat))
     bondz=np.zeros((max_n,nat))
    

     for i in range(nat):
       for j in range(i+1,nat):
         diff=coords[i]-coords[j]
         for v in range(3):
            if np.fabs(diff[v])>=dim05[v]:
                diff[v]=diff[v]-np.sign(diff[v])*dim[v]
                
         
         new_diff=np.dot(diff,lat)
         dd=np.linalg.norm(new_diff)
         bond_arr.append(dd)
     rcut=float(list(sorted(set(bond_arr)))[1]+list(sorted(set(bond_arr)))[2])/2.0
     #print 'rcut',rcut
     #print list(set(sorted(bond_arr)))
     bond_arr=[]
 
     for i in range(nat):
       for j in range(i+1,nat):
         diff=coords[i]-coords[j]
         for v in range(3):
            if np.fabs(diff[v])>=dim05[v]:
                diff[v]=diff[v]-np.sign(diff[v])*dim[v]
                
         
         new_diff=np.dot(diff,lat)
         dd=np.linalg.norm(new_diff)
         bond_arr.append(dd)
         if dd<rcut and dd>=0.1:
 
           nn_index=nn[i] #index of the neighbor
           nn[i]=nn[i]+1
           dist[nn_index][i]=dd #nn_index counter id
           nn_id[nn_index][i]=j #exact id
           bondx[nn_index,i]=new_diff[0]
           bondy[nn_index,i]=new_diff[1]
           bondz[nn_index,i]=new_diff[2]
           #print bond
           nn_index1=nn[j] #index of the neighbor
           nn[j]=nn[j]+1
           dist[nn_index1][j]=dd #nn_index counter id
           nn_id[nn_index1][j]=i #exact id
           bondx[nn_index1,j]=-new_diff[0]
           bondy[nn_index1,j]=-new_diff[1]
           bondz[nn_index1,j]=-new_diff[2]
     for i in range(nat):
      for in1 in range(nn[i]):
        j1=nn_id[in1][i]
        for in2 in range(in1+1,nn[i]):
          j2=nn_id[in2][i]
          nm=dist[in1][i]*dist[in2][i]
          if nm!=0:       
           cos=float(bondx[in1][i]*bondx[in2][i]+bondy[in1][i]*bondy[in2][i]+bondz[in1][i]*bondz[in2][i])/float(nm)
           if cos<=-1.0:
            cos=cos+0.000001
           if cos>=1.0:
            cos=cos-0.000001
           deg=math.degrees(math.acos(cos))
           deg_arr.append(deg)
          else:
            znm=znm+1
     print 'znm',znm 
     print 'deg',set(sorted(deg_arr))

     print 'dist',set(sorted(bond_arr))
     print 'nn',set(nn) 
     dist_hist, dist_bins = np.histogram(bond_arr,bins=np.arange(0, 10.0, .1), density=False)
     shell_vol = 4.0 / 3.0 * math.pi * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
     number_density = float(nat) /float(abs(np.dot(np.cross(lat[0], lat[1]), lat[2])))
     dist_hist = dist_hist / shell_vol/number_density#/len(nn)


     ang_hist, ang_bins = np.histogram(deg_arr,bins=np.arange(0, 180.0, 1), density=False)
     shell_vol = 4.0 / 3.0 * math.pi * (np.power(ang_bins[1:], 3) - np.power(ang_bins[:-1], 3))
     number_density = float(nat) /float(abs(np.dot(np.cross(lat[0], lat[1]), lat[2])))
     ang_hist = ang_hist / shell_vol/number_density#/len(nn)
     print len(dist_hist[:-1]) ,len(dist_bins)
     plt.plot(dist_bins[:-1],dist_hist)
     plt.xlabel('Distance (Ang)')
     plt.ylabel('Distribution')
     plt.savefig('rdf.png')
     plt.close()
     plt.plot(ang_bins[:-1],ang_hist)
     #plt.bar(ang_bins[:-1],ang_hist)
     plt.xlabel('Angle (degree)')
     plt.ylabel('Distribution')
     plt.savefig('angle.png')
     plt.close()


"""
p=Pos() #.read_poscar('POSCAR')
print p
p.verlet([4,4,4])
#p.write_poscar(fname='pp')
#print p.supercell([2,2,2]).frac_coords
#print p
#print p.frac_coords,p.elements,p.comment
import sys
sys.exit()
"""
