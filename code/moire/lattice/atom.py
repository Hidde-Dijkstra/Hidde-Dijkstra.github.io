import numpy as np
from dataclasses import dataclass
from typing import List

def rot_mat(θ):
    return np.array([[np.cos(θ), -np.sin(θ), 0], [np.sin(θ), np.cos(θ), 0], [0, 0, 1]])

class _Atom:
    
    lst = []
    
    def __init__(self, position, name='', draw_radius=0, weight=None, 
                 color=None, pseudo=dict(standard=[]), 
                 projections=[], scaling=1, **kwargs):
        self.id = str(len(_Atom.lst))
        _Atom.lst.append(self)
        self.position = np.array(position) / scaling
        self.name = name
        self.draw_radius = draw_radius / scaling
        self.pseudo = pseudo
        self.color = color
        self.weight = weight
        self.projections = projections
        self.NN = {}
        self.key = kwargs.get('key', name)
        
    def __eq__(self, other):
        return np.linalg.norm(self.position-other.position) < 10e-5
        
    def copy(self, position=None, scaling=1):
        if np.any(position == None):
            position = self.position
        return _Atom(position, name=self.name, draw_radius=self.draw_radius,
                     pseudo=self.pseudo, color=self.color, 
                     weight=self.weight, projections=self.projections,
                     key=self.key, scaling=scaling)
        
    def check_NN(self, other, R, lattice_vectors, 
                 metric=lambda x0, x1: np.linalg.norm(x0-x1)):
        if other.id not in self.NN:
            self.NN[other.id] = NN(atom=other, vectors=[], draw=True)
            other.NN[self.id] = NN(atom=self, vectors=[], draw=False)
        else:
            return
        for lat_vec in lattice_vectors:
            if metric(self.position, other.position+lat_vec) < R:
                self.NN[other.id].vectors.append(lat_vec)
                other.NN[self.id].vectors.append(-lat_vec)
                
@dataclass
class NN:
    atom: _Atom
    vectors: List[np.ndarray]
    draw: bool