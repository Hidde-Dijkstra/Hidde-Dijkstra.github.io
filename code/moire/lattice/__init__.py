import numpy as np
import plotly.graph_objects as go
import drawSvg as draw
from .atom import _Atom
from .plot import _DrawLattice


class Lattice:
    
    a_1 = np.array([1, 0, 0])
    a_2 = np.array([-1/2, np.sqrt(3)/2, 0])
    
    def __init__(self, a):
        self.a = a
        self.unit_cell = []
        self.atom_types = []

    def set_atom(self, position, atom_dic):
        if len(position) == 2:
            position = [*position, 0]
        self._add_atom(_Atom(position, **atom_dic))
        
    def _add_atom(self, atom, scale=True):
        if scale:
            atom = atom.copy(scaling=self.a)
        atom.NN = {}
        self.unit_cell.append(atom)
        if atom.name not in self.atom_types:
            self.atom_types.append(atom.name)
        
    def gen_NN(self, R, metric=lambda x0, x1: np.linalg.norm(x0-x1)):
        Δn = int(R/np.linalg.norm(self.a_1)) + 1
        n_range = range(-Δn, Δn+1)
        lattice_vectors = [self.a_1*i+self.a_2*j for i in n_range for j in n_range]
        for atom_1 in self.unit_cell:
            for atom_2 in self.unit_cell:
                atom_1.check_NN(atom_2, R, lattice_vectors, metric=metric)

    def plot(self, W, H=None, D=None, **kwargs):
        lat_plot = _DrawLattice(self, W, H, D)
        return lat_plot.plot(**kwargs)

class Supercell(Lattice):
    
    def __init__(self, lattice, m, n):
        scaling = np.sqrt(m**2+m*n+n**2)
        super().__init__(lattice.a*scaling)
        self.θ = np.arctan(np.sqrt(3)*(m+n)/(n-m))
        r = max(m, n)
        vec_range = [(self.a_1*i+self.a_2*j) for i in range(-2*r, 2*r) for j in range(-2*r, 2*r)]
        for atom in lattice.unit_cell:
            for vec in vec_range:
                atom_position = self._rot_mat(self.θ)@(atom.position+vec)
                if self._in_supercell(atom_position/scaling):
                    self._add_atom(atom.copy(atom_position, scaling=scaling), scale=False)

    def _rot_mat(self, θ):
        return np.array([[np.cos(θ), -np.sin(θ), 0], [np.sin(θ), np.cos(θ), 0], [0, 0, 1]])

    def _in_supercell(self, atom_position, tol=10**-5):
        M = np.linalg.inv(np.array([self.a_1[:2], self.a_2[:2]]).T)
        λ, μ = M @ atom_position[:2]
        in_parallellogram = (-tol < λ < 1-tol) and (-tol < μ < 1-tol) 
        return in_parallellogram

class Bilayer(Lattice):
    
    def __init__(self, layer_1, layer_2):
        super().__init__(layer_1.a)
        self.θ = (np.pi-np.abs(layer_1.θ-layer_2.θ))/np.pi*180
        for layer in [layer_1, layer_2]:
            for atom in layer.unit_cell:
                self._add_atom(atom.copy(atom.position), scale=False)