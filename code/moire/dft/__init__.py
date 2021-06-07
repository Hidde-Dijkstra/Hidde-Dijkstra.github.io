from .decorators import check
from moire.bandstructure import BandStructure
import numpy as np

from .wannier import Wannier
from .w90 import W90
from .qe import QE
from .projwfc import PROJWFC

class DFT:

    a_1 = np.array([1, 0, 0])
    a_2 = np.array([-1/2, np.sqrt(3)/2, 0])
    a_3 = np.array([0, 0, 1])

    def __init__(self, prefix, work_dir, data_dir='', Δz=32):
        self.prefix = prefix
        self.work_directory = work_dir
        self.data_directory = data_dir
        self.qe = False
        self.w90 = False
        self.Δz = Δz

    def prepare_qe(self, qe_dic, lattice):
        self.lattice = lattice
        self.atoms = {}
        for atom in lattice.unit_cell:
            if atom.name not in self.atoms:
                self.atoms[atom.name] = atom.__dict__
                self.atoms[atom.name]['loc'] = [atom.position*lattice.a]
            else:
                self.atoms[atom.name]['loc'].append(atom.position*lattice.a)
        self.qe = True
        self.QE = QE(self, qe_dic, lattice)
        return self.QE

    @check('qe')
    def prepare_projwfc(self, projwfc_dir=''):
        self.PROJWFC = PROJWFC(self, projwfc_dir)
        return self.PROJWFC
    
    @check('qe')
    def prepare_w90(self, w90_dic):
        self.w90 = True
        self.W90 = W90(self, w90_dic)
        return self.W90

    @check('w90')
    def wannier_orbitals(self, orbital_dir=''):
        self.Wannier = Wannier(self, orbital_dir)
        return self.Wannier

    def set_k_path(self, k_list, k_tags, N_k):
        self.k_list = k_list
        self.k_tags = k_tags
        self.N_k = N_k
        for i in range(len(k_list)):
            self.k_list[i] = np.array(k_list[i])
        bs = self.gen_bs()
        self.spacing = bs.spacing
        self.k_spacing = bs.k_spacing
        self.k_path = bs.k_path  
    
    def set_k_grid(self, N):
        self.N_grid = N 
        self.k_grid = [np.array([i/N, j/N, 0]) for i in range(N) for j in range(N)] 

    def gen_bs(self):
        bs = BandStructure()
        bs.set_k_path(self.k_list, self.k_tags, self.N_k)
        return bs 