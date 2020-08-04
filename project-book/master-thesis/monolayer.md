---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.9'
    jupytext_version: 1.5.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Hoi boi

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
from numpy.linalg import matrix_power as mat_pow
import matplotlib.pyplot as plt
from datetime import datetime
plt.rcParams['figure.figsize'] = [15, 10]
```

```{code-cell} ipython3
:tags: [remove-output, hide-input]

class LatVec:
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
        
    def rot(self, θ):
        R_θ = np.array([[np.cos(θ), -np.sin(θ)], [np.sin(θ), np.cos(θ)]])
        return R_θ @ self.vec
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
        self.NN_array = np.zeros((self.N_atoms, 6, 2), dtype=int)
        for i in range(self.N_atoms):
            self.NN_array[i, :] = [self.find_NN(i, h_vec) for h_vec in self.hop_list]
            
    def find_NN(self, i, h_vec):
        for m, lat_vec in enumerate((self.v_1 & self.v_2) + [LatVec(0, 0)]):
            if self.grid[i]+h_vec+lat_vec in self.grid:
                return self.grid.index(self.grid[i]+h_vec+lat_vec), m
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

bla bla

```{glue:figure} my-fig
---
figclass: margin-caption
name: figfig
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer accumsan tellus nec pulvinar lacinia. Pellentesque sed tristique augue. Aenean in congue lorem. Duis commodo auctor est quis pharetra. Aenean interdum at ligula id blandit. Vestibulum feugiat dignissim leo, eu maximus ex convallis sed. Sed dolor lectus, viverra non nulla et, hendrerit vulputate nisi. Etiam nec orci laoreet, facilisis nunc id, semper ex.
```
bla bla

```{figure} ../logo.png
---
figclass: margin-caption
name: figfig2
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer accumsan tellus nec pulvinar lacinia. Pellentesque sed tristique augue. Aenean in congue lorem. Duis commodo auctor est quis pharetra. Aenean interdum at ligula id blandit. Vestibulum feugiat dignissim leo, eu maximus ex convallis sed. Sed dolor lectus, viverra non nulla et, hendrerit vulputate nisi. Etiam nec orci laoreet, facilisis nunc id, semper ex.
```

+++

Dit was de functie:

$$
f(x) = x^2
$$

```{code-cell} ipython3
from bokeh.plotting import figure, show, output_notebook

# prepare some data
x = [1, 2, 3, 4, 5]
y = [6, 7, 2, 4, 5]


# create a new plot with a title and axis labels
p = figure(title="simple line example", x_axis_label='x', y_axis_label='y')

# add a line renderer with legend and line thickness
p.line(x, y, legend_label="Temp.", line_width=2)
output_notebook()
# show the results
show(p)

glue("my-fig2", p, display=True)
```

```{glue:figure} my-fig2
---
figclass: margin-caption
name: figfig2
---
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer accumsan tellus nec pulvinar lacinia. Pellentesque sed tristique augue. Aenean in congue lorem. Duis commodo auctor est quis pharetra. Aenean interdum at ligula id blandit. Vestibulum feugiat dignissim leo, eu maximus ex convallis sed. Sed dolor lectus, viverra non nulla et, hendrerit vulputate nisi. Etiam nec orci laoreet, facilisis nunc id, semper ex.
```

```{code-cell} ipython3

```
