from jarvis.core.lattice import Lattice
from jarvis.core.atoms import Atoms
import numpy as np
from numpy import linalg as LA
from collections import OrderedDict
import pprint

class Kpoints(object):
    def __init__(self,kpoints=[[1,1,1]],labels=[],kpoint_mode='automatic',header='Gamma'):
       self._kpoints = kpoints
       self._labels = labels
       self._kpoint_mode = kpoint_mode
       self._header = header
    def automatic_length_mesh(self,lattice_mat=[],length=20, header='Gamma'):
        inv_lat = Lattice(lattice_mat=lattice_mat).inv_lattice()      
        b1 = LA.norm(np.array(inv_lat[0])) 
        b2 = LA.norm(np.array(inv_lat[1])) 
        b3 = LA.norm(np.array(inv_lat[2])) 
        n1=int(max(1, length * b1 + 0.5))
        n2=int(max(1, length * b2 + 0.5))
        n3=int(max(1, length * b3 + 0.5))
        return Kpoints(kpoints=[[n1,n2,n3]],header=header,kpoint_mode='automatic')
    
    def write_file(self,filename=''):
        if self._kpoint_mode=='automatic':
           f=open(filename,'w')
           f.write('Automatic kpoint scheme\n')    
           f.write('0\n')
           line=str(self._header)+'\n'
           f.write(line)
           kp = self._kpoints
           line=str(kp[0][0])+str(' ')+str(kp[0][1])+str(' ')+str(kp[0][2])+str('\n')
           f.write(line)
           f.close()
          
       
    def to_dict(self):
         d = OrderedDict()
         d['kpoints']=self._kpoints
         d['labels']=self._labels 
         d['kpoint_mode']=self._kpoint_mode
         d['header']=self._header
         return d
 
    def __repr__(self,indent=4):
         return pprint.pformat(self.to_dict(), indent=indent) 



class HighSymmetry3DKpointFactory(object):
    def __init__(self,kpoints=[],path=[]):
        self._kpoints=kpoints
        self._path=path
    def cubic(self):
        """
        Cubic HighSymmKPath
        :return: Dict
        """
        self.name = "CUB"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "R", "X"], ["M", "R"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def fcc(self):
        """
        Fcc HighSymmKPath
        :return: Dict
        """
        self.name = "FCC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'K': np.array([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'U': np.array([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
                   'W': np.array([0.5, 1.0 / 4.0, 3.0 / 4.0]),
                   'X': np.array([0.5, 0.0, 0.5])}
        path = [["\\Gamma", "X", "W", "K",
                 "\\Gamma", "L", "U", "W", "L", "K"], ["U", "X"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def bcc(self):
        """
        Bcc HighSymmKPath
        :return: Dict
        """
        self.name = "BCC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'H': np.array([0.5, -0.5, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "H", "N", "\\Gamma", "P", "H"], ["P", "N"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def tet(self):
        """
        Tetragonal HighSymmKPath
        :return: Dict
        """
        self.name = "TET"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "R", "A", "Z"], ["X", "R"],
                ["M", "A"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def bctet1(self, c, a):
        """
        BCT1 HighSymmKPath
        :return: Dict
        """
        self.name = "BCT1"
        eta = (1 + c ** 2 / a ** 2) / 4.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'M': np.array([-0.5, 0.5, 0.5]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'Z': np.array([eta, eta, -eta]),
                   'Z_1': np.array([-eta, 1 - eta, eta])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "P", "N", "Z_1", "M"],
                ["X", "P"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def bctet2(self, c, a):
        """
        BCT2 HighSymmKPath
        :return: Dict
        """
        self.name = "BCT2"
        eta = (1 + a ** 2 / c ** 2) / 4.0
        zeta = a ** 2 / (2 * c ** 2)
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   '\\Sigma': np.array([-eta, eta, eta]),
                   '\\Sigma_1': np.array([eta, 1 - eta, -eta]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'Y': np.array([-zeta, zeta, 0.5]),
                   'Y_1': np.array([0.5, 0.5, -zeta]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        path = [["\\Gamma", "X", "Y", "\\Sigma", "\\Gamma", "Z",
                 "\\Sigma_1", "N", "P", "Y_1", "Z"], ["X", "P"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def orc(self):
        """
        Orthorhombic HighSymmKPath
        :return: Dict
        """
        self.name = "ORC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'S': np.array([0.5, 0.5, 0.0]),
                   'T': np.array([0.0, 0.5, 0.5]),
                   'U': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "S", "Y", "\\Gamma",
                 "Z", "U", "R", "T", "Z"], ["Y", "T"], ["U", "X"], ["S", "R"]]

        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)
    def orcf1(self, a, b, c):
        """
        Orthorhombic f1 HighSymmKPath
        :return: Dict
        """
        self.name = "ORCF1"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4

        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5 + zeta, zeta]),
                   'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'T': np.array([1, 0.5, 0.5]),
                   'X': np.array([0.0, eta, eta]),
                   'X_1': np.array([1, 1 - eta, 1 - eta]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                ["T", "X_1"], ["X", "A", "Z"], ["L", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def orcf2(self, a, b, c):
        """
        Orthorhombic f2 HighSymmKPath
        :return: Dict
        """
        self.name = "ORCF2"
        phi = (1 + c ** 2 / b ** 2 - c ** 2 / a ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        delta = (1 + b ** 2 / a ** 2 - b ** 2 / c ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'C': np.array([0.5, 0.5 - eta, 1 - eta]),
                   'C_1': np.array([0.5, 0.5 + eta, eta]),
                   'D': np.array([0.5 - delta, 0.5, 1 - delta]),
                   'D_1': np.array([0.5 + delta, 0.5, delta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'H': np.array([1 - phi, 0.5 - phi, 0.5]),
                   'H_1': np.array([phi, 0.5 + phi, 0.5]),
                   'X': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "C", "D", "X", "\\Gamma",
                 "Z", "D_1", "H", "C"], ["C_1", "Z"], ["X", "H_1"], ["H", "Y"],
                ["L", "\\Gamma"]]

        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)
    def orcf3(self, a, b, c):
        """
        Orthorhombic f3 HighSymmKPath
        :return: Dict
        """
        self.name = "ORCF3"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5 + zeta, zeta]),
                   'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'T': np.array([1, 0.5, 0.5]),
                   'X': np.array([0.0, eta, eta]),
                   'X_1': np.array([1, 1 - eta, 1 - eta]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                ["X", "A", "Z"], ["L", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def orci(self, a, b, c):
        """
        Orthorhombic I HighSymmKPath
        :return: Dict
        """
        self.name = "ORCI"
        zeta = (1 + a ** 2 / c ** 2) / 4
        eta = (1 + b ** 2 / c ** 2) / 4
        delta = (b ** 2 - a ** 2) / (4 * c ** 2)
        mu = (a ** 2 + b ** 2) / (4 * c ** 2)
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([-mu, mu, 0.5 - delta]),
                   'L_1': np.array([mu, -mu, 0.5 + delta]),
                   'L_2': np.array([0.5 - delta, 0.5 + delta, -mu]),
                   'R': np.array([0.0, 0.5, 0.0]),
                   'S': np.array([0.5, 0.0, 0.0]),
                   'T': np.array([0.0, 0.0, 0.5]),
                   'W': np.array([0.25, 0.25, 0.25]),
                   'X': np.array([-zeta, zeta, zeta]),
                   'X_1': np.array([zeta, 1 - zeta, -zeta]),
                   'Y': np.array([eta, -eta, eta]),
                   'Y_1': np.array([1 - eta, eta, -eta]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        path = [["\\Gamma", "X", "L", "T", "W", "R", "X_1", "Z",
                 "\\Gamma", "Y", "S", "W"], ["L_1", "Y"], ["Y_1", "Z"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def orcc(self, a, b, c):
        """
        Orthorhombic C HighSymmKPath
        :return: Dict
        """
        self.name = "ORCC"
        zeta = (1 + a ** 2 / b ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([zeta, zeta, 0.5]),
                   'A_1': np.array([-zeta, 1 - zeta, 0.5]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   'T': np.array([-0.5, 0.5, 0.5]),
                   'X': np.array([zeta, zeta, 0.0]),
                   'X_1': np.array([-zeta, 1 - zeta, 0.0]),
                   'Y': np.array([-0.5, 0.5, 0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "S", "R", "A", "Z",
                 "\\Gamma", "Y", "X_1", "A_1", "T", "Y"], ["Z", "T"]]

        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)
    def hex(self):
        """
        Hexagonal HighSymmKPath
        :return: Dict
        """
        self.name = "HEX"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'H': np.array([1.0 / 3.0, 1.0 / 3.0, 0.5]),
                   'K': np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]),
                   'L': np.array([0.5, 0.0, 0.5]),
                   'M': np.array([0.5, 0.0, 0.0])}
        path = [["\\Gamma", "M", "K", "\\Gamma", "A", "L", "H", "A"], ["L", "M"],
                ["K", "H"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def rhl1(self, alpha):
        """
        Rhombohedral 1 HighSymmKPath
        :return: Dict
        """
        self.name = "RHL1"
        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'B': np.array([eta, 0.5, 1.0 - eta]),
                   'B_1': np.array([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0]),
                   'L_1': np.array([0.0, 0.0, -0.5]),
                   'P': np.array([eta, nu, nu]),
                   'P_1': np.array([1.0 - nu, 1.0 - nu, 1.0 - eta]),
                   'P_2': np.array([nu, nu, eta - 1.0]),
                   'Q': np.array([1.0 - nu, nu, 0.0]),
                   'X': np.array([nu, 0.0, -nu]),
                   'Z': np.array([0.5, 0.5, 0.5])}
        path = [["\\Gamma", "L", "B_1"], ["B", "Z", "\\Gamma", "X"],
                ["Q", "F", "P_1", "Z"], ["L", "P"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def rhl2(self, alpha):
        """
        Rhombohedral 2 HighSymmKPath
        :return: Dict
        """
        self.name = "RHL2"
        eta = 1 / (2 * tan(alpha / 2.0) ** 2)
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([0.5, -0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0]),
                   'P': np.array([1 - nu, -nu, 1 - nu]),
                   'P_1': np.array([nu, nu - 1.0, nu - 1.0]),
                   'Q': np.array([eta, eta, eta]),
                   'Q_1': np.array([1.0 - eta, -eta, -eta]),
                   'Z': np.array([0.5, -0.5, 0.5])}
        path = [["\\Gamma", "P", "Z", "Q", "\\Gamma",
                 "F", "P_1", "Q_1", "L", "Z"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mcl(self, b, c, beta):
        """
        Monoclinic 1 HighSymmKPath
        :return: Dict
        """
        self.name = "MCL"
        eta = (1 - b * cos(beta) / c) / (2 * sin(beta) ** 2)
        nu = 0.5 - eta * c * cos(beta) / b
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5, 0.0]),
                   'C': np.array([0.0, 0.5, 0.5]),
                   'D': np.array([0.5, 0.0, 0.5]),
                   'D_1': np.array([0.5, 0.5, -0.5]),
                   'E': np.array([0.5, 0.5, 0.5]),
                   'H': np.array([0.0, eta, 1.0 - nu]),
                   'H_1': np.array([0.0, 1.0 - eta, nu]),
                   'H_2': np.array([0.0, eta, -nu]),
                   'M': np.array([0.5, eta, 1.0 - nu]),
                   'M_1': np.array([0.5, 1 - eta, nu]),
                   'M_2': np.array([0.5, 1 - eta, nu]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'Y': np.array([0.0, 0.0, 0.5]),
                   'Y_1': np.array([0.0, 0.0, -0.5]),
                   'Z': np.array([0.5, 0.0, 0.0])}
        path = [["\\Gamma", "Y", "H", "C", "E", "M_1", "A", "X", "H_1"],
                ["M", "D", "Z"], ["Y", "D"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mclc1(self, a, b, c, alpha):
        """
        Monoclinic C1 HighSymmKPath
        :return: Dict
        """
        self.name = "MCLC1"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                   'F_1': np.array([zeta, zeta, eta]),
                   'F_2': np.array([-zeta, -zeta, 1 - eta]),
                   # 'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                   'I': np.array([phi, 1 - phi, 0.5]),
                   'I_1': np.array([1 - phi, phi - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([1 - psi, psi - 1, 0.0]),
                   'X_1': np.array([psi, 1 - psi, 0.0]),
                   'X_2': np.array([psi - 1, -psi, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'Y_1': np.array([-0.5, -0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                ["Y", "X_1"], ["X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mclc2(self, a, b, c, alpha):
        """
        Monoclinic C2 HighSymmKPath
        :return: Dict
        """
        self.name = "MCLC2"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                   'F_1': np.array([zeta, zeta, eta]),
                   'F_2': np.array([-zeta, -zeta, 1 - eta]),
                   'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                   'I': np.array([phi, 1 - phi, 0.5]),
                   'I_1': np.array([1 - phi, phi - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([1 - psi, psi - 1, 0.0]),
                   'X_1': np.array([psi, 1 - psi, 0.0]),
                   'X_2': np.array([psi - 1, -psi, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'Y_1': np.array([-0.5, -0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                ["N", "\\Gamma", "M"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mclc3(self, a, b, c, alpha):
        """
        Monoclinic C3 HighSymmKPath
        :return: Dict
        """
        self.name = "MCLC3"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                   'F_1': np.array([phi, phi - 1, psi]),
                   'F_2': np.array([1 - phi, -phi, 1 - psi]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([0.5, -0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "H", "Z", "I", "F_1"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mclc4(self, a, b, c, alpha):
        """
        Monoclinic C4 HighSymmKPath
        :return: Dict
        """
        self.name = "MCLC4"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                   'F_1': np.array([phi, phi - 1, psi]),
                   'F_2': np.array([1 - phi, -phi, 1 - psi]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([0.5, -0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "H", "Z", "I"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def mclc5(self, a, b, c, alpha):
        """
        Monoclinic C5 HighSymmKPath
        :return: Dict
        """
        self.name = "MCLC5"
        zeta = (b ** 2 / a ** 2 + (1 - b * cos(alpha) / c)
                / sin(alpha) ** 2) / 4
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        mu = eta / 2 + b ** 2 / (4 * a ** 2) - b * c * cos(alpha) / (2 * a ** 2)
        nu = 2 * mu - zeta
        rho = 1 - zeta * a ** 2 / b ** 2
        omega = (4 * nu - 1 - b ** 2 * sin(alpha) ** 2 / a ** 2) * c / (2 * b * cos(alpha))
        delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([nu, nu, omega]),
                   'F_1': np.array([1 - nu, 1 - nu, 1 - omega]),
                   'F_2': np.array([nu, nu - 1, omega]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([rho, 1 - rho, 0.5]),
                   'I_1': np.array([1 - rho, rho - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "H", "F_1"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def tria(self):
        """
        Trigonal a HighSymmKPath
        :return: Dict
        """
        self.name = "TRI1a"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.0]),
                   'M': np.array([0.0, 0.5, 0.5]),
                   'N': np.array([0.5, 0.0, 0.5]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)

    def trib(self):
        """
        Trigonal b HighSymmKPath
        :return: Dict
        """
        self.name = "TRI1b"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, -0.5, 0.0]),
                   'M': np.array([0.0, 0.0, 0.5]),
                   'N': np.array([-0.5, -0.5, 0.5]),
                   'R': np.array([0.0, -0.5, 0.5]),
                   'X': np.array([0.0, -0.5, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.0]),
                   'Z': np.array([-0.5, 0.0, 0.5])}
        path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
        return HighSymmetry3DKpointFactory(kpoints=kpoints,path=path)
          
if __name__=='__main__':
   box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
   coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
   elements = ["Si", "Si"]
   Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
   lattice_mat=Si.lattice_mat
   kp=Kpoints().automatic_length_mesh(lattice_mat=lattice_mat)
   print (kp.__repr__(0)) 
   kp.write_file('KPOINTS')
