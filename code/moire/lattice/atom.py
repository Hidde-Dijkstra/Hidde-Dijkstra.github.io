import numpy as np
import plotly.graph_objects as go

def rot_mat(θ):
    return np.array([[np.cos(θ), -np.sin(θ), 0], [np.sin(θ), np.cos(θ), 0], [0, 0, 1]])

class _Atom:
    
    lst = []
    
    def __init__(
            self, 
            position,
            name='', 
            draw_radius=0, 
            pseudo_potential=None, 
            color=None, 
            weight=None,
            scaling=1, 
            projections=[],
            key=''
    ):
        self.id = str(len(_Atom.lst))
        _Atom.lst.append(self)
        self.position = np.array(position)
        self.name = name
        self.draw_radius = draw_radius
        self.pseudo_potential = pseudo_potential
        self.color = color
        self.weight = weight
        self.projections = projections
        self.NN = {}
        self.key = key
        
    def copy(self, position=None, scaling=1):
        if np.any(position == None):
            position = self.position
        return _Atom(
            position / scaling, 
            name=self.name, 
            draw_radius=self.draw_radius / scaling, 
            pseudo_potential=self.pseudo_potential, 
            color=self.color, 
            weight=self.weight, 
            projections=self.projections,
            key=self.key
        )
        
    def check_NN(self, other, R, lattice_vectors, metric=lambda x0, x1: np.linalg.norm(x0-x1)):
        if other.id not in self.NN:
            self.NN[other.id] = dict(atom=other, vecs=[], draw=True)
            other.NN[self.id] = dict(atom=self, vecs=[], draw=False)
        else:
            return
        for lat_vec in lattice_vectors:
            if metric(self.position, other.position+lat_vec) < R:
                self.NN[other.id]['vecs'].append(lat_vec)
                other.NN[self.id]['vecs'].append(-lat_vec)