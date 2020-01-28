from jarvis.io.vasp.outputs import Vasprun

Ry_to_ev = 13.6056980659
class WriteInputs(object):
      

      def __init__(self,vasprun_path='', energy=None, struct=None, intrans=None):


           self.energy=energy
           self.struct=struct
           self.intrans=intrans
           self.vasprun_path=vasprun_path
           self.vrun=Vasprun(filename=vasprun_path)

      def write_energy(self,filename='boltztrap.energyso',trim=0.1):
          kpoints=self.vrun.kpoints._kpoints
          eigs_up,eigs_dn=self.vrun.eigenvalues
          ef=self.vrun.efermi     
          target = 2*int(len(eigs_dn[0])*(1-trim))#+1
          print ('target',target)
          f=open(filename,'w')
          line=str('system \n')+str(len(kpoints))+'\n'
          f.write(line)
          for i,j,k in zip(kpoints,eigs_up,eigs_dn):
            count = 0
            line=(' '.join(map(str,i)))+str(' ')+str(target)+'\n'
            #f.write(line)
            f.write("%12.8f %12.8f %12.8f %d\n"%(i[0],i[1],i[2],target))
            for m,n in zip(j,k):
             count = count+2
             if count<=target:
              en_up = round((m[0]-ef)/float(Ry_to_ev),8)
              en_dn = round((n[0]-ef)/float(Ry_to_ev),8)
              f.write("%18.8f\n" % (en_up))
              f.write("%18.8f\n" % (en_dn))
          f.close() 
if __name__ == "__main__" :
 from jarvis.io.vasp.outputs import Vasprun
 vrun=Vasprun('/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/vasprun.xml')
 print (vrun.final_energy)
 inp=WriteInputs(vasprun_path='/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/vasprun.xml')
 inp.write_energy()
