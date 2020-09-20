---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.12'
    jupytext_version: 1.5.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import drawSvg as draw
from myst_nb import glue
import physdraw as phd
```

```{code-cell} ipython3
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
    
    def __eq__(self, other):
        return (self.i==other.i) & (self.j==other.j)
    
    def __and__(self, other):
        # A simple way to get all hopping vectors from two unit vectors.
        return [self, other.neg(), self.neg()+other.neg(), self.neg(), other, self+other]
    
    def __mul__(self, other):
        if type(other) == LatVec:
            return np.dot(self.vec, other.vec)
        else:
            return np.dot(self.vec, other)
        
    def __rmul__(self, other):
        return self * other
    
    def neg(self):
        return LatVec(-self.i, -self.j)
    
    def vectorize(self):
        if self.reciprocal:
            return self.scale*(self.i*self.b_1 + self.j*self.b_2)
        else:
            return self.scale*(self.i*self.a_1 + self.j*self.a_2)
    
    def plot(self, container, atom_radius=0.2, atom_color='darkblue', θ=0, atom="", bonds=False, displace=0, **kwargs):
        origin = self.vec + displace
        a_list = [self.a_1, self.a_2, self.a_1-self.a_2]      
        if θ != 0:
            origin = rot_mat(θ) @ self.vec
            a_list = [rot_mat(θ) @ a for a in a_list]
        if bonds: 
            for a in a_list:
                container.append(draw.Line(*(origin-a/2), *(origin+a/2), **kwargs))
        gradient = draw.RadialGradient(*origin, atom_radius)
        gradient.addStop(0, 'white', 1)
        gradient.addStop(1, atom_color, 1)
        container.append(draw.Circle(*origin, atom_radius, fill=gradient, **kwargs))
        container.append(draw.Text(atom, atom_radius, *origin, text_anchor='middle', alignment_baseline="central"))
        return container
    
container = draw.Drawing(4, 2.3, origin='center', displayInline=False)
for (i, j) in [(i, j) for i in range(-3, 3) for j in range(-2, 2)]:
    container = LatVec(i, j).plot(container, stroke='black', atom_color='red', bonds=True, atom=str(i)+','+str(j), stroke_width=0.015)
    container = LatVec(i, j).plot(container, stroke='black', atom_color='blue', atom="Se2", stroke_width=0.01, 
                                  displace=[1/2, 1/3], atom_radius=0.15)
container.setRenderSize(1000)
glue('fig:lattice', container, display=False)
```

```{code-cell} ipython3


container = draw.Drawing(1000, 350, origin="center")
d_z2(container)
d_xy(container, translate=(300, 0))
d_xy(container, translate=(-300, 0))
glue('fig:d_z2-d_xy', container, display=False)
```

```{code-cell} ipython3
var_dic = {
    "t_1": 0.034, 
    "t_2": 0.263, 
    "t_3": -0.207, 
    "t_12": 0.329, 
    "t_13": 0.486, 
    "t_23": 0.457, 
    "ε_1": 2.179, 
    "ε_3": 0.942, 
    "λ_SOC": 0.228
}
t_1, t_2, t_3, t_12, t_13, t_23, ε_1, ε_3, λ_SOC = var_dic.values()
for var_key in var_dic:
    glue("var:"+var_key, var_dic[var_key], display=False)
```

```{code-cell} ipython3
class Supercell:
    hop_list = LatVec(1, 0) & LatVec(0, -1)
    
    def __init__(self, m, n):
        self.v_1 = LatVec(m, n)    
        self.v_2 = LatVec(n+m, -m)
        self.w_1 = LatVec(m, n+m, True, 1/(m**2+m*n+n**2))
        self.w_2 = LatVec(n, -m, True, 1/(m**2+m*n+n**2))
        r = max(m, n)
        self.grid = [LatVec(i, j) for i in range(0, 3*r) for j in range(-r, r+1) if self.in_supercell(i, j)]
        self.N_atoms = len(self.grid)
        self.Δθ = np.arctan((n-m)/(n+m)/np.sqrt(3))
        self.construct_NN_array()
         
    def in_supercell(self, i, j, tol=10**-5):
        M = np.linalg.inv(np.array([self.v_1.vec, self.v_2.vec]).T)
        λ, μ = M @ LatVec(i, j).vec
        in_parellogram = (tol < λ < 1-tol) and (tol < μ < 1-tol) 
        return in_parellogram or (i, j) == (0, 0)
    
    def construct_NN_array(self):
        self.NN_array = []
        for i in range(self.N_atoms):
            self.NN_array.append([self.find_NN(i, h_vec) for h_vec in self.hop_list])
            
    def find_NN(self, i, h_vec):
        for lat_vec in (self.v_1 & self.v_2) + [LatVec(0, 0)]:
            if self.grid[i]+h_vec+lat_vec in self.grid:
                return {'NN': self.grid.index(self.grid[i]+h_vec+lat_vec),
                        'LatVec': lat_vec
                       }
        raise Exception('No NN found for '+str(i)+' '+str(h_vec))
                
    def interlayer_hopping_array(self, supercell, tol=10**-5):
        if self.N_atoms != supercell.N_atoms:
            raise Exception('Supercells have a different number of atoms')
        if np.abs(self.Δθ + supercell.Δθ) > tol:
            raise Exception('Unequal twist angles')
        z_hopping = np.zeros((self.N_atoms, self.N_atoms, 2))
        for i in range(self.N_atoms):
            vec_i = supercell.grid[i].rot(-supercell.Δθ)    
            for j in range(self.N_atoms):
                min_ΔR = 10**6
                for lat_vec in (self.v_1 & self.v_2) + [LatVec(0, 0)]:
                    vec_j_trial = (self.grid[j]+lat_vec).rot(-self.Δθ)
                    if np.linalg.norm(vec_i-vec_j_trial) < min_ΔR:
                        min_ΔR = np.linalg.norm(vec_i-vec_j_trial)
                        vec_j = vec_j_trial
                z_hopping[i, j] = vec_i - vec_j
        return z_hopping
                        
    def plot(self, W, H, rotate=True, atom_color='blue'):
        if rotate == True:
            supercell = phd.Lattice(self.v_1.vec, self.v_2.vec, W, H, θ=self.Δθ)
        else:
            supercell = phd.Lattice(self.v_1.vec, self.v_2.vec, W, H)        
        for atom in self.grid:
            supercell.add_atom(phd.LatticeAtom(atom.vec, atom_color=atom_color))
        return supercell.draw_lattice()
        
    def plot_supercell(self, grid_points=None, *, lat_vec=LatVec(0, 0), rotate=True, **kwargs):
        if grid_points == None:
            grid_points = range(self.N_atoms)
        grid_array = np.array([(self.grid[i]+lat_vec).rot(-self.Δθ*rotate) for i in grid_points])
        plt.scatter(grid_array[:, 0], grid_array[:, 1], **kwargs)
```

```{code-cell} ipython3
n, m = 6, 7
layer_1 = Supercell(n, m)
layer_2 = Supercell(m, n)
z_hopping = layer_1.interlayer_hopping_array(layer_2)
layer_1.Δθ*2*180/np.pi

j_1 = 2
j_2 = 44
for lat_vec in [layer_1.v_1+layer_1.v_2, layer_1.v_1,  layer_1.v_2, LatVec(0, 0)]:
    layer_1.plot_supercell(lat_vec=lat_vec, c='green', s=10)
    layer_1.plot_supercell(list(layer_1.NN_array[j_1, :]), lat_vec=lat_vec, c='blue', marker='d')
    layer_1.plot_supercell([j_1], lat_vec=lat_vec, c='blue', s=70, marker='D')
    layer_1.plot_supercell([j_2], lat_vec=lat_vec, c='blue', marker='D', s=70)
for lat_vec in [layer_2.v_1+layer_2.v_2, layer_2.v_1,  layer_2.v_2, LatVec(0, 0)]:
    layer_2.plot_supercell(lat_vec=lat_vec, c='orange', s=10)
    layer_2.plot_supercell([i for i in range(layer_2.N_atoms) if np.linalg.norm(z_hopping[i, j_2]) < 1], lat_vec=lat_vec, c='red', marker='d')

plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')
```

```{code-cell} ipython3

```

```{code-cell} ipython3
0.1 == 0.1
```

```{code-cell} ipython3

```
