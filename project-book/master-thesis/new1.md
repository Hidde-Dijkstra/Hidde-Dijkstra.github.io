---
jupytext:
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

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import drawSvg as draw
from myst_nb import glue
import physdraw as phd
import moire_supercell as moire
from numpy.linalg import matrix_power as mat_pow
```

```{code-cell} ipython3
class WSe2:
    E_mat = np.array([[0.034, 0.329, 0.486],
                      [-0.329, 0.263, 0.457],
                      [0.486, -0.457, -0.207]])
    E_levels = np.diag([2.179, 2.179, 0.943])
    E_0 = E_levels[0, 0] - 3/2*(E_mat[1, 1]+E_mat[0, 0]) - 3*np.sqrt(3)*E_mat[0, 1]
    E_1 = E_levels[2, 2] - 3*E_mat[2, 2]
    λ_SOC = 0.228
    d = 2.4
    λ = 0.3
    
    def __init__(self, m, n, R_max=1):
        self.layer_1 = moire.Supercell(m, n)
        self.layer_2 = moire.Supercell(n, m)
        self.z_hopping = self.layer_1.interlayer_hopping_array(self.layer_2)
        R = np.array([[-1/2, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -1/2, 0], [0, 0, 1]])
        C_2 = np.diag([1, -1, 1])
        self.E_list = [mat_pow(R, i)@mat_pow(C_2, i)@self.E_mat@(mat_pow(R, i)@mat_pow(C_2, i)).T for i in range(6)]
        z_hopping = self.layer_1.interlayer_hopping_array(self.layer_2)
        self.z_hop = []
        for i in range(self.layer_1.N_atoms):
            self.z_hop += [[(j, z_hopping[i, j]) for j in range(self.layer_2.N_atoms) if np.linalg.norm(z_hopping[i, j])<R_max]]
    
        
    def H_WSe2_mono(self, k, σ, layer, rotate=False):
        if rotate:
            k = np.array([[np.cos(layer.Δθ), -np.sin(layer.Δθ)], [np.sin(layer.Δθ), np.cos(layer.Δθ)]]) @ k
        H = np.zeros((3*layer.N_atoms, 3*layer.N_atoms), dtype=complex)
        for i in range(layer.N_atoms):
            H[3*i:3*(i+1), 3*i:3*(i+1)] = self.E_levels + self.λ_SOC*np.array([[0, 1j*σ, 0], [-1j*σ, 0, 0], [0, 0, 0]])
            for j in range(6):
                m = layer.NN_array[i][j]['NN']
                H[3*i:3*(i+1), 3*m:3*(m+1)] = self.E_list[j] * np.exp(1j*np.dot(layer.NN_array[i][j]['LatVec'].vec, k))
        return H
    
    def plot_mono_bands(self, σ, k_list, N_points, layer=None, k_tags=None, **kwargs):
        if layer == None:
            layer = self.layer_1
        spacing, k_path_unit = self.k_path(k_list, N_points)
        N_points = len(k_path_unit)
        k_path_super = self.k_path([layer.reduce_k_point(k) for k in k_list], N_points, spacing)[1]
        #k_path_super = k_path_unit
        ε_supercell = np.zeros((N_points, 3*layer.N_atoms))
        ε_unitcell = np.zeros((N_points, 3))
        for i in range(N_points):
            my_H = (sum([self.E_list[j] * np.exp(1j*np.dot(layer.hop_list[j].vec, k_path_unit[i])) for j in range(6)]) 
                    + self.E_levels + self.λ_SOC*np.array([[0, 1j*σ, 0], [-1j*σ, 0, 0], [0, 0, 0]]))
            ε_unitcell[i, :] = np.linalg.eigh(my_H)[0]
            ε_supercell[i, :] = np.linalg.eigh(self.H_WSe2_mono(k_path_super[i], σ, layer))[0]
        for x in [sum(spacing[:i]) for i in range(len(spacing)+1)]:
            plt.axvline(x=x, linewidth=0.5, color='black', ls='--')
        plt.axhline(y=0, linewidth=0.5, color='black', ls='--')
        plt.axhline(y=self.E_1-self.E_0, linewidth=0.5, color='black', ls='--')
        plt.plot(range(N_points), ε_supercell-self.E_0, c='b')
        plt.plot(range(N_points), ε_unitcell-self.E_0, c='r')
        plt.ylabel('Energy (eV)')
        if k_tags != None:
            plt.xticks([sum(spacing[:i]) for i in range(len(spacing)+1)], k_tags)
            
    def plot_bilayer_bands(self, σ, t_0, k_list, N_points, k_tags=None, print_time=False, **kwargs):
        spacing, k_path_super = self.k_path([self.layer_1.reduce_k_point(k) for k in k_list], N_points)
        N_points = len(k_path_super)
        ε_supercell = np.zeros((N_points, 6*self.layer_1.N_atoms))
        for i in range(N_points):
            startTime = datetime.now()
            ε_supercell[i, :] = np.linalg.eigh(self.H_WSe2_bilayer(k_path_super[i], σ, t_0))[0]
            if print_time:
                print(datetime.now()-startTime)
        for x in [sum(spacing[:i]) for i in range(len(spacing)+1)]:
            plt.axvline(x=x, linewidth=0.5, color='black', ls='--')
        plt.axhline(y=0, linewidth=0.5, color='black', ls='--')
        plt.axhline(y=self.E_1-self.E_0, linewidth=0.5, color='black', ls='--')
        plt.plot(range(N_points), ε_supercell-self.E_0, c='b')
        plt.ylabel('Energy (eV)')
        if k_tags != None:
            plt.xticks([sum(spacing[:i]) for i in range(len(spacing)+1)], k_tags)
        
    
    def k_path(self, k_list, N_points, spacing=None):
        if spacing == None:
            k_norms = [np.linalg.norm(k_list[i+1]-k_list[i]) for i in range(len(k_list)-1)]
            spacing = [int(N_points*k_norms[i]/sum(k_norms)) for i in range(len(k_norms))]
        k_path = []
        for i in range(len(spacing)):
            k_path += [k_list[i] + (k_list[i+1]-k_list[i])*j/spacing[i] for j in range(spacing[i])]
        k_path += [k_list[-1]]
        return spacing, k_path
    
    def H_WSe2_bilayer(self, k, σ, t_0):
        H = np.zeros((6*self.layer_1.N_atoms, 6*self.layer_1.N_atoms), dtype=complex)
        H[0:3*self.layer_1.N_atoms, 0:3*self.layer_1.N_atoms] = self.H_WSe2_mono(k, σ, self.layer_1, rotate=True)
        H[3*self.layer_1.N_atoms:, 3*self.layer_1.N_atoms:] = self.H_WSe2_mono(k, σ, self.layer_2, rotate=True)
        for i in range(self.layer_1.N_atoms):
            for j, hop_vec in self.z_hop[i]:
                t = t_0 * np.exp(-(np.sqrt(np.dot(hop_vec, hop_vec)+self.d**2)-self.d)/self.λ)
                H[3*i+2+3*self.layer_1.N_atoms, 3*j+2] = t * np.exp(1j*np.dot(k, hop_vec))
                H[3*j+2, 3*i+2+3*self.layer_1.N_atoms] = t * np.exp(-1j*np.dot(k, hop_vec))
        return H
        
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

```

```{code-cell} ipython3
0.1 == 0.1
```

```{code-cell} ipython3
container = draw.Drawing(1000, 350, origin="center")
d_z2(container)
d_xy(container, translate=(300, 0))
d_xy(container, translate=(-300, 0))
glue('fig:d_z2-d_xy', container, display=False)
```
