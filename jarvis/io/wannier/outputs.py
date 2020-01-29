import sys
import matplotlib
matplotlib.use('Agg') #fixes display issues?
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math
import time
import copy as copy
from jarvis.core.kpoints import Kpoints3D
from jarvis.io.vasp.inputs import Poscar

class WannierHam(object):
    def __init__(self,filename='wannier90_hr.dat',nwan=None,nr=None, sym_r=None, H_int=None, H_val=None, H=None, HR=None):
         self.filename=filename
         self.nr=nr
         self.nrwan=nwan
         self.sym_r=sym_r
         self.H_int=H_int
         self.H_val=H_val
         self.H=H
         self.HR=HR
  
         if self.nr is None:
             self.read_ham()

    def read_ham(self):
        f=open(self.filename,'r')
        lines = f.read().splitlines()
        allines=f.readlines()
        f.close()
        self.nwan=int(lines[1])
        self.nr=int(lines[2])
                
        lines_r = int(math.floor(self.nr / 15) + 1)      

        #print ('self.nwan,self.nr',self.nwan,self.nr,lines_r)
        self.sym_r = np.zeros(self.nr, dtype=float)
        #load sym ops
        for i in range(lines_r):

            num = 3+i
            start = i*15
            end = (i+1)*15
            if end > self.nr:
                end = self.nr
            self.sym_r[start:end] = ([float(j) for j in lines[num].split()])
        #print (self.sym_r)

        tot = self.nwan**2 * self.nr
        self.H_int = np.zeros((tot,5),dtype=int)
        self.H_val = np.zeros(tot,dtype=complex)


        c=0
        rnum = 0
        for i in range(lines_r+3,lines_r+3+tot):

            rnum = (c)//self.nwan**2
            #print ('lines[i].split()[0:5]',c,self.nwan**2)#,self.sym_r[rnum])
            self.H_int[c,:] = [int(j) for j in lines[i].split()[0:5]]
            #print ('lines[i][5]',lines[i])
            self.H_val[c] = float(lines[i].split()[5]) + 1j * float(lines[i].split()[6]) #/ float(self.sym_r[rnum])
            c+=1

        #print (self.H_int[0,:])
        #print (self.H_val[0])

        #print (self.H_int[-1,:])
        #print (self.H_val[-1])

        #print ('loaded ', self.filename)
        #print ('nwan: ', self.nwan)
        #print ('nr:   ', self.nr)
        #print ()


        #reshape

        nx1 = np.min(self.H_int[:,0])
        nx2 = np.max(self.H_int[:,0])
        ny1 = np.min(self.H_int[:,1])
        ny2 = np.max(self.H_int[:,1])
        nz1 = np.min(self.H_int[:,2])
        nz2 = np.max(self.H_int[:,2])

        self.ind = [[nx1,nx2],[ny1,ny2],[nz1,nz2]]

        ix = nx2-nx1+1
        iy = ny2-ny1+1
        iz = nz2-nz1+1

        print ('H size', ix,iy,iz,self.nwan,self.nwan)

        self.H = np.zeros((ix,iy,iz,self.nwan,self.nwan),dtype=complex)

        self.ind_dict = {}
        for i in range(self.H_val.shape[0]):
            ind = self.get_ind(self.H_int[i,0:3])

            self.ind_dict[(ind[0],ind[1],ind[2])] = i

            nw1 = self.H_int[i,3]
            nw2 = self.H_int[i,4]

            self.H[ind[0], ind[1], ind[2], nw1-1,nw2-1] = self.H_val[i]

        print ('done reshaping1')
        nr = ix*iy*iz
        self.R = np.zeros((nr,3),dtype=float)
        self.HR = np.zeros((nr,self.nwan**2),dtype=complex)

        c=0
        for x in range(nx1, nx2+1):
            for y in range(ny1, ny2+1):
                for z in range(nz1, nz2+1):
                    #                    ind = self.ind_dict[(x,y,z)]
                    ind = self.get_ind([x,y,z])
                    self.R[c,:] = [x,y,z]
                    self.HR[c,:] = np.reshape(self.H[ind[0],ind[1],ind[2], :,:], self.nwan**2)
                    c+=1

        if c != nr:
            ValueError('c is not equal to r', c, nr)

        return self.R, self.H, self.HR


    def get_ind(self,nxyz):

        return [nxyz[0] - self.ind[0][0], nxyz[1] - self.ind[1][0], nxyz[2] - self.ind[2][0]]


    def solve_ham(self,k=[0,0,0], proj=None):


        print ('solve', self.nwan, self.R.shape, self.HR.shape)

        nr = self.R.shape[0]
        print ('nr==',nr,self.nr)
        hk = np.zeros((self.nwan,self.nwan),dtype=complex)


        kmat = np.tile(k, (nr,1))

        exp_ikr = np.exp(-1.0j*2*np.pi* np.sum(kmat*self.R, 1))

        temp = np.zeros(self.nwan**2, dtype=complex)
        for i in range(nr):
            temp += exp_ikr[i]*self.HR[i,:]

        hk = np.reshape(temp, (self.nwan, self.nwan))

        hk = (hk + hk.T.conj())/2.0

        val, vect = np.linalg.eigh(hk)

        if proj is not None:
            p = np.real(np.sum(vect[proj,:]*np.conj(vect[proj, :]), 0))
        else:
            p = np.ones(val.shape)

#        print 'proj', np.sum(np.sum(p))

        return val.real, vect, p

    def band_structure_eigs(self,kpath=None,proj=None,efermi=0.0):  
       eigs=[]
       for i in kpath:
          print (i)
          val, vect,p = self.solve_ham(k=i,proj=proj)
          eigs.append(val-efermi)
       return np.array(eigs)


    def get_bandstructure_plot(self,atoms=None,efermi=0.0,filename='bs.png',yrange=[-4,4]):
         kpoints,labels=Kpoints3D().interpolated_points(atoms) #_path

         #print ('kpath',kpoints)

         eigs=self.band_structure_eigs(kpath=kpoints,efermi=efermi).T

         for i,ii in enumerate(eigs):
            plt.plot(ii,color='b')

         plt.ylim(yrange)

         plt.savefig(filename)

         plt.close()

                  
        
class Wannier90wout(object):
        def __init__(self, wout_path='wannier90.wout'):
           self.wout = wout_path

        def give_wannier_centers(self):
         f=open(self.wout,'r')
         lines=f.read().splitlines()
         f.close()
         final=False
         wan_cnts=[]
         for ii,i in enumerate(lines):
          if 'Final State' in i:
            final=True
          if final:
            if 'WF centre and spread' in i:
             tmp=[float(j) for j in i.split('(')[1].split(')')[0].split(',')]
             wan_cnts.append(tmp)
         return wan_cnts

if __name__ == "__main__":
  hr = '/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/wannier90_hr.dat'
  wout = '/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/wannier90.wout'
  centers =  Wannier90wout(wout_path=wout).give_wannier_centers()
  #print (centers)

  #val,vect,p = WannierHam(filename=hr).solve_ham()
  #print ('val=',val)
  #print ('vect=',vect)
  #print ('p=',p)
 
  p=Poscar.from_file('/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/POSCAR').atoms
  print (p)
  WannierHam(filename=hr).get_bandstructure_plot(atoms=p)

  #R,H, HR = WannierHam(filename=hr).read_ham()
  #print ('R',R)
  #print ()
  #print ()
  #print ('HR',HR)
