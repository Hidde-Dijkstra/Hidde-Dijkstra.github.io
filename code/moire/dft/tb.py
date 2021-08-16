from .decorators import change_directory, save
from . import io_processing as io
from dataclasses import dataclass
from moire.plot.heatmap import gen_heatmap
from collections import defaultdict
from moire.plot.figure import figure
import numpy as np
from scipy.interpolate import griddata
import plotly.graph_objs as go

class TB:
    
    def __init__(self, DFT):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.a = DFT.QE.qe_dic['&SYSTEM']['a']
        self.a_vec = [self.a*DFT.a_1, self.a*DFT.a_2, self.a*DFT.a_3]
        self.keys = []
      
    @change_directory('data_directory')  
    def orbitals(self, use_wan=False):
        wan_centers = iter(np.load(f'{self.prefix}/wan_centers.npy'))
        return [Orbital(next(wan_centers) if use_wan else self.a*atom.position, 
                        projection, i, spin=1-2*σ)
                for i, atom in enumerate(self.DFT.lattice.unit_cell)
                for projection in atom.projections
                for σ in range(1+self.DFT.W90.w90_dic['spinors'])]

    @change_directory('data_directory')
    @save
    @change_directory('work_directory')       
    def extract(self):
        hopping_file = io.read_file(f'{self.prefix}_hr.dat')
        hop_list = [io.gen_lst(lst, '   ', float) 
                    for lst in io.gen_lst(hopping_file, '\n')
                    if len(io.gen_lst(lst, '   ')) == 7]
        return dict(tb_hopping=hop_list)
    
    def H(self, R_max, print_N=True):
        N = self.DFT.W90.w90_dic['num_wann']
        hop_dic = defaultdict(lambda: np.zeros((N, N), dtype=complex))
        N_hop = 0
        orbitals = self.orbitals()
        for n_vec, i, j, t in self.get_hoppings():
            ΔR = np.dot(n_vec, self.a_vec) - (orbitals[i]-orbitals[j])*self.a
            if np.dot(ΔR[:2], ΔR[:2]) <= R_max**2:
                N_hop += 1
                hop_dic[n_vec][i, j] += t
        def H(k):
            return sum([hop_mat*np.exp(2*np.pi*1j*np.dot(n_vec, k))
                        for n_vec, hop_mat in hop_dic.items()])
        if print_N:
            print(N_hop, 'hoppings found')
        return H
    
    def plot_potential(self, spin=False, atom_list=None, use_wan=False, 
                       **kwargs):
        hop_dic = defaultdict(list)
        orbitals = self.orbitals(use_wan)
        def on_site(i, j, n_vec):
            return i == j and not np.any(n_vec)
        for n_vec, i, j, t in self.get_hoppings(on_site):
            key = f'{orbitals[i]} ' + (str(orbitals[i].spin) if spin else '')
            if atom_list == None or orbitals[i].i in atom_list:
                hop_dic[key].append(np.array([*orbitals[i].position[:2], 
                                              np.abs(t)]))
        slider = dict(args=dict(positions=list(hop_dic.values())), 
                      tags=list(hop_dic.keys()))
        return gen_heatmap(slider=slider, a=self.a, **kwargs)
    
    def plot_interlayer(self, radial=True, use_wan=False, spinful=False, **kwargs):
        coupling = defaultdict(list)
        orbitals = self.orbitals(use_wan)
        def interlayer(i, j, _):
            return i < len(orbitals)/2 and j >= len(orbitals)/2
        for n_vec, i, j, t in self.get_hoppings(interlayer):
            if not spinful or (i-j)%2 != 0:
                key, value = self.ΔR(n_vec, i, j, t, orbitals, radial=radial)
                coupling[key].append(value)
        kwargs.update(dict(slider=dict(args=dict(coupling_lst=list(coupling.values())),
                                       tags=list(coupling.keys())),
                           showlegend=False))      
        if radial:
            fig = self._plot_radial_interlayer(**kwargs)
            fig.update_yaxes(title='θ')
        else:
            fig = self._plot_interlayer(**kwargs)
            fig.update_yaxes(scaleanchor="x", title='Å', scaleratio=1)
        fig.update_xaxes(title='Å', constrain='domain')
        return fig
        
    @figure
    def _plot_radial_interlayer(self, coupling_lst=None, R_max=2, N=100, 
                                **kwargs):
        R, θ, t = np.transpose([(σ*R, θ+2*np.pi*i, t) 
                                for i in [-1, 0, 1] 
                                for σ in [-1, 1]
                                for R, θ, t in coupling_lst])
        Ri, θi = np.meshgrid(np.linspace(0, R_max, N), 
                             np.linspace(-np.pi, np.pi, N))
        return [go.Heatmap(x=np.linspace(0, R_max, N), 
                           y=np.linspace(-np.pi, np.pi, N),
                           z=griddata(np.transpose([R, θ]), np.abs(t), (Ri, θi)))]
        
    @figure
    def _plot_interlayer(self, coupling_lst=None, x_max=2, y_max=1, N=100,
                         **kwargs):
        x, y, t = np.transpose(coupling_lst)
        xi, yi = np.meshgrid(np.linspace(-x_max, x_max, N), 
                             np.linspace(-y_max, y_max, N))
        return [go.Heatmap(x=np.linspace(-x_max, x_max, N), 
                           y=np.linspace(-y_max, y_max, N),
                           z=griddata(np.transpose([x, y]), np.abs(t), (xi, yi)))]
    
    @change_directory('data_directory')
    def get_hoppings(self, check=lambda i, j, n_vec: True):
        return [(tuple(hop[:3]), int(hop[3])-1, int(hop[4])-1, hop[5]+1j*hop[6])
                for hop in np.load(f'{self.prefix}/tb_hopping.npy')
                if check(int(hop[3])-1, int(hop[4])-1, tuple(hop[:3]))]
    
    def ΔR(self, n_vec, i, j, t, orbitals, radial=True, spin=False):
        ΔR = np.dot(n_vec, self.a_vec) - (orbitals[i]-orbitals[j])
        key = f'{orbitals[i]}-{orbitals[j]}'
        if f'{orbitals[j]}-{orbitals[i]}' in self.keys and i != j:
            key = f'{orbitals[j]}-{orbitals[i]}'
            σ = -1
        else:
            key = f'{orbitals[i]}-{orbitals[j]}'
            σ = 1
            if key not in self.keys:
                self.keys.append(key)                
        if radial:
            θ = np.angle(ΔR[0]+1j*ΔR[1])
            if σ == -1:
                θ = θ - np.pi if θ > 0 else θ + np.pi
            value = (np.abs(ΔR[0]+1j*ΔR[1]), θ, np.abs(t))
        else:
            value = (σ*ΔR[0], σ*ΔR[1], np.abs(t))
            if value[2] > 0.02:
                print(value, orbitals[i].position-orbitals[j].position, orbitals[i], orbitals[j])
        return key, value

@dataclass
class Orbital:
    position: np.ndarray
    orbital: str
    i: int 
    spin: int
    
    def __sub__(self, other):
        return self.position - other.position
    
    def __str__(self):
        return self.orbital