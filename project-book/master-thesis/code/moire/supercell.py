# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Supercell

from moire import diagram_gen as dg
import numpy as np


class LatVec:
    # We identify each atom as the two integers i and j which connect it to the origin. 
    # Using a pythonic object we can define how two of these vectors interact.
    a_1 = np.array([1, 0])
    a_2 = np.array([1/2, np.sqrt(3)/2])
    b_1 = 2*np.pi * np.array([1, -1/np.sqrt(3)])
    b_2 = 2*np.pi * np.array([0, 2/np.sqrt(3)])
    
    
    def __init__(self, i, j, reciprocal=False, scale=1):
        self.i = i
        self.j = j
        self.scale = scale
        self.reciprocal = reciprocal
        self.vec = self.vectorize()
          
    def __add__(self, other):
        return LatVec(self.i+other.i, self.j+other.j)
    
    def __sub__(self, other):
        return LatVec(self.i-other.i, self.j-other.j)
    
    def __eq__(self, other):
        return (self.i==other.i) & (self.j==other.j)
    
    def __and__(self, other):
        # A simple way to get all hopping vectors from two unit vectors.
        return [self, other.neg(), self.neg()+other.neg(), self.neg(), other, self+other]
    
    def __mul__(self, other):
        if type(other) == LatVec:
            return np.dot(self.vec, other.vec)
        elif type(other) == int:
            return LatVec(self.i*other, self.j*other)
        else:
            return np.dot(self.vec, other)
        
    def __rmul__(self, other):
        return self * other
    
    def neg(self):
        return LatVec(-self.i, -self.j)
    
    def rot(self, θ):
         return np.array([[np.cos(θ), -np.sin(θ)], [np.sin(θ), np.cos(θ)]]) @ self.vec
    
    def vectorize(self):
        if self.reciprocal:
            return self.scale*(self.i*self.b_1 + self.j*self.b_2)
        else:
            return self.scale*(self.i*self.a_1 + self.j*self.a_2)


class Supercell:
    supercell_hop = [(0, 0), (1, 0), (0, -1), (-1, -1), (-1, 0), (0, 1), (1, 1)]
    atom_hop = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

    def __init__(self, m, n):
        self.v = [LatVec(m, n), LatVec(n+m, -m)]
        self.w_1 = LatVec(m, n+m, True, 1/(m**2+m*n+n**2))
        self.w_2 = LatVec(n, -m, True, 1/(m**2+m*n+n**2))
        r = max(m, n)
        self.grid = [LatVec(i, j) for i in range(0, 3*r) for j in range(-r, r+1) if self.in_supercell(i, j)]
        self.N_atoms = len(self.grid)
        self.Δθ = np.arctan((n-m)/(n+m)/np.sqrt(3))
        self.intralayer_NN = []
        for i in range(self.N_atoms):
            self.intralayer_NN.append([self.find_NN(i, LatVec(*h_vec)) for h_vec in self.atom_hop])
        
    def in_supercell(self, i, j, tol=10**-5):
        M = np.linalg.inv(np.array([v.vec for v in self.v]).T)
        λ, μ = M @ LatVec(i, j).vec
        in_parellogram = (tol < λ < 1-tol) and (tol < μ < 1-tol) 
        return in_parellogram or (i, j) == (0, 0)
        
    def find_NN(self, i, h_vec):
        for hop in self.supercell_hop:
            lat_vec = dg.Index(*hop) * self.v
            if self.grid[i]+h_vec-lat_vec in self.grid:
                return [self.grid.index(self.grid[i]+h_vec-lat_vec), *hop]
        raise Exception('No NN found for '+str(i)+' '+str(h_vec))
                                        
    def lattice(self, W, H, rotate=True, NN_intralayer=False, atom_color='blue', atom_radius=0.2):
        supercell = dg.Lattice(*[v.vec for v in self.v], W, H, θ=self.Δθ*rotate)
        for atom in self.grid:
            supercell.add_atom(dg.LatticeAtom(atom.vec, atom_color=atom_color, atom_radius=atom_radius))
        if NN_intralayer:
            for i, atom in enumerate(supercell.unit_cell):
                for NN in self.intralayer_NN[i]:
                        supercell.NN(atom, supercell.unit_cell[NN[0]], [NN[1:]])                                    
        return supercell
    
    def reduce_k_point(self, k):
        α, β = np.linalg.inv([[self.w_1*self.w_1, self.w_1*self.w_2], [self.w_1*self.w_2, self.w_2*self.w_2]]) @ np.array([np.dot(k, self.w_1.vec), np.dot(k, self.w_2.vec)]).T
        return np.modf(α)[0]*self.w_1.vec + np.modf(β)[0]*self.w_2.vec


class Bilayer:
    hop_list = [(0, 0), (1, 0), (0, -1), (-1, -1), (-1, 0), (0, 1), (1, 1)]
    
    
    def __init__(self, m, n, tol=10**-5):
        self.layer_1 = Supercell(m, n)
        self.layer_2 = Supercell(n, m)
        if self.layer_1.N_atoms != self.layer_2.N_atoms:
            raise Exception('Supercells have a different number of atoms')
        self.N =  self.layer_1.N_atoms
        if np.abs(self.layer_1.Δθ + self.layer_2.Δθ) > tol:
            raise Exception('Unequal twist angles')
        self.v = [dg.rot_mat(-self.layer_1.Δθ) @ v.vec for v in self.layer_1.v]
        self.interlayer_NN = self.find_interlayer_NN()
               
    def find_interlayer_NN(self, max_distance=1, tol=10**-5):
        interlayer_NN = []
        for i in range(self.N):
            interlayer_NN.append([])
            vec_i = self.layer_1.grid[i].rot(-self.layer_1.Δθ)    
            for j in range(self.N):
                min_ΔR = 10**6
                for hop in self.hop_list:
                    lat_vec = dg.Index(*hop) * self.layer_2.v
                    vec_j_trial = (self.layer_2.grid[j]+lat_vec).rot(-self.layer_2.Δθ)
                    if np.linalg.norm(vec_i-vec_j_trial) < min_ΔR:
                        min_ΔR = np.linalg.norm(vec_i-vec_j_trial)
                        vec_j = vec_j_trial
                        my_hop = hop
                if np.linalg.norm(vec_i - vec_j) < max_distance+tol:
                    interlayer_NN[i].append([j, *my_hop])
        return interlayer_NN
    
    def lattice(self, W, H, atom_color=['blue', 'red'], atom_radius=0.2, NN_interlayer=False, NN_intralayer=False, add_Se2=False,
               Se2_color=['blueviolet', 'coral']):
        bilayer = dg.Lattice(*self.v, W, H)
        Se2_list = []
        for i, layer in enumerate([self.layer_1, self.layer_2]):
            rot = dg.rot_mat(-layer.Δθ)
            for atom in layer.grid:
                bilayer.add_atom(dg.LatticeAtom(rot@atom.vec, atom_color=atom_color[i], atom_radius=atom_radius))
                if add_Se2:
                    pos = atom.vec + (1/2, (1-2*i)/(2*np.sqrt(3)))
                    Se2_list.append(dg.LatticeAtom(rot@pos, atom_color=Se2_color[i], atom_radius=atom_radius/2))
        if NN_interlayer:
            for i, atom in enumerate(bilayer.unit_cell[:self.N]):
                for NN in self.interlayer_NN[i]:
                    bilayer.NN(atom, bilayer.unit_cell[self.N+NN[0]], [NN[1:]])                    
                    bilayer.NN(bilayer.unit_cell[self.N+NN[0]], atom, [[-NN[1], -NN[2]]])
        if NN_intralayer:
            for j, layer in enumerate([self.layer_1, self.layer_2]):
                for i, atom in enumerate(bilayer.unit_cell[j*self.N:(j+1)*self.N]):
                    for NN in layer.intralayer_NN[i]:
                        bilayer.NN(atom, bilayer.unit_cell[self.N*j+NN[0]], [NN[1:]])
        if add_Se2:
            for Se2 in Se2_list:
                bilayer.add_atom(Se2)
        return bilayer


