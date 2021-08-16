from .decorators import check
from moire.bandstructure import Path
import numpy as np
from .wannier import Wannier
from .w90 import W90
from .qe import QE
from .relax import RELAX
from .projwfc import PROJWFC
from .tb import TB

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

    def prepare_qe(self, qe_dic, lattice, middle=False):
        self.middle = middle
        self.lattice = lattice
        self.qe = True
        self.QE = QE(self, qe_dic)
        return self.QE
    
    @check('qe')
    def prepare_relax(self, relax_dir=''):
        self.RELAX = RELAX(self, relax_dir)
        return self.RELAX

    @check('qe')
    def prepare_projwfc(self, projwfc_dir=''):
        self.PROJWFC = PROJWFC(self, projwfc_dir)
        return self.PROJWFC
    
    @check('qe')
    def prepare_w90(self, w90_dic, N_super=3):
        self.w90 = True
        self.W90 = W90(self, w90_dic, N_super=N_super)
        return self.W90
    
    @check('w90')
    def tight_binding(self):
        self.TB = TB(self)
        return self.TB

    @check('w90')
    def wannier_orbitals(self, orbital_dir=''):
        self.Wannier = Wannier(self, orbital_dir)
        return self.Wannier

    def set_k_path(self, k_list, tags, N):
        self.path = Path(k_list, tags, N)
    
    def set_k_grid(self, N):
        self.N_grid = N 
        self.k_grid = [np.array([i/N, j/N, 0]) for i in range(N) for j in range(N)] 