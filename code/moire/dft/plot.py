import numpy as np
from .decorators import change_directory, check
import plotly.graph_objects as go
from moire import lattice as dg
from moire.plot import figure, BandStructure
from scipy.interpolate import RectBivariateSpline

@change_directory('data_directory')
def plot_bands(self, ticks=None):
    if self.qe:
        self.band_plot.add_bands(np.load('qe_bands.npy'), 'blue', 'qe bands')
    if self.w90:
        self.band_plot.fig.update_yaxes(
            range=[self.w90_dic['dis_win_min'], self.w90_dic['dis_win_max']]
        )
        self.band_plot.add_bands(
            np.load('w90_bands.npy'), 
            'green', 
            'w90 bands', 
            k_arr=np.load('w90_band_ticks.npy')
        )

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

@change_directory('data_directory')
def plot_wfc(self, k_list=[], k_tags=[], N_k=[], **kwargs):
    bs = BandStructure(cut=kwargs.get('cut'))
    bs.set_k_path(k_list, k_tags, N_k)
    bs.add_bands(bands=np.load(self.prefix+'/projwfc/bands.npy'))
    return _plot_wfc(self=self, fig=bs.fig, **kwargs)

@change_directory('data_directory')
@figure
def _plot_wfc(
        self, ε_range=[-np.inf, np.inf], orbital={}, 
        elements={}, size=15,**kwargs):
    bands = np.load(self.prefix+'/projwfc/bands.npy')
    states = sort_states(np.load(self.prefix+'/projwfc/states.npy'))
    occupations = np.load(self.prefix+'/projwfc/occupations.npy')
    occupations = occupations / np.expand_dims(np.sum(occupations, axis=2), axis=2)

    wfc = {}
    for key in ['atoms', 'l', 'm_l', 'm_s']:
        if key in orbital:
            wfc[key] = list_transform(orbital[key])
        elif key in elements:
            wfc[key] = list_transform(elements[key])
        else:
            wfc[key] = []

    orbital_states = check_states(states, wfc)
    orbital_density = np.sum(occupations[:, :, orbital_states], axis=2)
    traces = []

    for i, band in enumerate(bands.T):
        if np.logical_and(ε_range[0]<band, band<ε_range[1]).any():
            marker_dic = dict(
                size=size*orbital_density[:, i],
                line=dict(width=1),
            )
            x = np.arange(len(band))
            traces.append(go.Scattergl(
                x = x, 
                y = band,
                opacity=0.3,
                text = orbital_density[:, i],
                mode = 'markers', 
                marker = marker_dic
            ))
    return traces
    
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

def sort_states(states):
    atom, l, j, m_j = states.T
    m_l = m_j - np.sign((j-l)*m_j) / 2
    m_s = np.sign(j-l)*np.sign(m_j) / 2
    return np.array([atom, l, m_l, m_s])

def check_states(states, wfc):
    valid_states = []
    for i, state in enumerate(states.T):
        valid_state = True
        for key, value in zip(['atoms', 'l', 'm_l', 'm_s'], state):
            valid_state &= (np.abs(value) in wfc[key]) or (wfc[key] == [])
        if valid_state:
            valid_states.append(i)
    return valid_states

def list_transform(variable):
    if type(variable) == list:
        return variable
    else:
        return [variable]