import xmltodict
import types
from pymatgen.io.vasp.outputs import Vasprun

class Vasprun(object):

      def __init__(self,filename='vasprun.xml',data= {}):
         self._filename = filename
         self._data = data
         self.ionic_steps = None
         self.electronic_steps = None
         if self._data=={}:
            self.xml_to_dict()


      def xml_to_dict(self):
          with open(self._filename) as fd:
              data = xmltodict.parse(fd.read())
              self._data = data
              self.ionic_steps = data['modeling']['calculation']
              if type(self.ionic_steps) is not list:
                    self.ionic_steps = [self.ionic_steps]
              
      @property
      def final_energy(self):
          return float(self.ionic_steps[-1]['scstep'][-1]['energy']['i'][11]['#text'])



      @property
      def num_atoms(self):
          return self._data['modeling']['atominfo']['atoms']

      @property
      def num_types(self):
          return int(self._data['modeling']['atominfo']['types'])

      @property
      def elements(self):
          elements = [self._data['modeling']['atominfo']['array'][0]['set']['rc'][i]['c'][0] for i in range(len(self._data['modeling']['atominfo']['array'][0]['set']['rc']))]
          if len(elements)!=self.num_atoms:
             ValueError ('Number of atoms is  not equal to number of elements')
          return elements
if __name__=='__main__':
    filename='/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/vasprun.xml'
    v=Vasprun(filename=filename)
    #print (v._filename, v.final_energy)
    print (v._filename,'elements', v.elements)
    #print ('pmg',Vasprun(v._filename).final_energy)
    filename='/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-BAND-bulk@mp_541837/vasprun.xml'
    v=Vasprun(filename=filename)
    #print (v._filename,v.final_energy)
    print (v._filename,v.elements)
    #print ('pmg',Vasprun(v._filename).final_energy)
    #print (v._data.keys())
    #print ()
    #print ()
    #print (v._data['modeling'].keys())    
    #print ()
    #print ()
    ##'generator', 'incar', 'kpoints', 'parameters', 'atominfo', 'structure', 'calculation'
    #print ('incar',v._data['modeling']['incar'])
    #print ()
    #print ()
    #print ('kpoints',v._data['modeling']['kpoints'])
    #print ()
    #print ()
    #print ('atominfo',v._data['modeling']['atominfo'])
    #print ()
    #print ()
    #print ('parameters',v._data['modeling']['parameters'])
    #print ()
    #print ()
    #print ('structure',v._data['modeling']['structure'])
    #print ()
    #print ()
    #print ('calculation',v._data['modeling']['calculation'].keys())
    #print ()
    #print ()
    ##calculation odict_keys(['scstep', 'structure', 'varray', 'energy', 'time', 'eigenvalues', 'separator', 'dos', 'projected'])
