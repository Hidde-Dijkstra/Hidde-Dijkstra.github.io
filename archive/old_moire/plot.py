import numpy as np
from .decorators import change_directory, check
import plotly.graph_objects as go
from moire import lattice as dg
from moire.plot import figure, BandStructure
from scipy.interpolate import RectBivariateSpline
from sympy.physics.quantum.cg import CG

@check('qe')
@change_directory('data_directory')
def plot_nscf(self):
    nscf_data = np.tile(np.load(self.prefix+'/qe_nscf.npy')[:64, :], (3, 3, 1))
    x = np.linspace(-1, 2, len(nscf_data))
    bands = np.zeros((len(self.band_plot.k_path), nscf_data.shape[2]))
    k_x, k_y, _ = np.transpose(self.band_plot.k_path)
    for i in range(nscf_data.shape[2]):
        f_interpolation = RectBivariateSpline(x, x, nscf_data[:, :, i])
        bands[:, i] = f_interpolation(k_x, k_y, grid=False)
    self.band_plot.add_bands(bands, 'blue', 'nscf')

@change_directory('data_directory')
def plot_relax(self):
    pass
    
def interpolate(data):
    x = np.linspace(-1, 2, 3*len(data))
    return RectBivariateSpline(x, x, np.tile(data, (3, 3, 1)))

def plot_tight_binding(
        self, 
        R_max, 
        name='tb', 
        color='red', 
        show_amount_NN=False
    ):
    self.band_plot.plot_from_H(
        self.create_H(R_max, show_amount_NN=show_amount_NN), 
            name, 
            color, 
            N_input=3
        )
            
def convert_grid(A):
    n_1, n_2, n_3 = A.shape
    grid = np.zeros((n_1+int(n_2/2)-1, int(n_2/2), n_3))
    for i in range(n_1):
        for j in range(int(n_2/2)):
            grid[int(n_2/2)-1+i-j, j, :] = A[i, 2*j, :]
    n = int(n_3/2)
    Δn = int(n_2/2)
    return grid[::2, :, n-Δn:n-Δn+n_2]




def untangle2(unsorted_bands, α=1):
    M, N = unsorted_bands.shape
    bands = np.zeros((M, N))
    bands[:2] = unsorted_bands[:2]
    δ_bands = bands[1] - bands[0]
    for m in range(2, M):
        δ_bands = α*(bands[m-1]-bands[m-2]) + (1-α)*δ_bands     
        bands_estimate = δ_bands + bands[m-1]
        bands = bands[:, np.argsort(bands_estimate)]
        δ_bands = δ_bands[np.argsort(bands_estimate)]
        bands[m] = unsorted_bands[m]
    return bands