import numpy as np
from .decorators import change_directory, check
import plotly.graph_objects as go
from moire import lattice as dg
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
def plot_wfc(
        self,
        ε_range, 
        orbital=None, 
        i_orbital=0, 
        size=15, 
        showlegend=True
    ):
    bands = np.load(self.prefix+'/projwfc/bands.npy')
    states = sort_states(np.load(self.prefix+'/projwfc/states.npy'))
    occupations = np.load(self.prefix+'/projwfc/occupations.npy')
    occupations /= np.sum(occupations, axis=2)

    orbital_states = orbital.states(states, i_orbital)
    orbital_density = np.sum(occupations[:, :, orbital_states], axis=2)

    traces = []

    for band in bands.T:
        if np.logical_and(ε_range[0]<band, band<ε_range[1]).any():
            marker_dic = dict(
                size=size*orbital_density,
                color=orbital.colors[i_orbital],
                line=dict(width=1),
            )
            x = np.arange(len(band))
            traces.append(go.Scattergl(
                x = x, 
                y = band,
                opacity=0.3,
                legendgroup = orbital.name[i_orbital], 
                text = orbital_density,
                name = orbital.name[i_orbital],
                mode = 'markers', 
                showlegend = showlegend, 
                marker = marker_dic
            ))
            showlegend = False
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
    atom, l, j, m_j = states
    m_l = m_j - np.sign((j-l)*m_j) / 2
    m_s = np.sign(j-l)*np.sign(m_j) / 2
    return np.array([atom, l, m_l, m_s])

class Orbital:

        def __init__(
            self, 
            names=[''], 
            colors=['black'],
            layer=False,
            spin=False, 
            layer_1=[], 
            layer_2=[]
        ):
            self.names = names
            self.colors = colors
            self.layers = [layer_1, layer_2]
            self.type_list = []
            if layer:
                pass
            self.spin = spin
            self.options = [[layer, spin] for spin in []]

        def check_state(states, i):
            for state in states:
                atom, l, m_l, m_s = state
                include = atom in []

        
def layout(func):
    def wrapper(**kwargs):
        

def legend(func):
    def wrapper(**kwargs):
        if 'legend' in kwargs:
            legend = kwargs['legend']
            key, args = list(legend['args'].items())[0]
            traces = []
            for i in range(len(args)):
                traces.append(func(**{key: args[i]}, **kwargs))
            return traces
        else:
            return func(**kwargs)
    return wrapper

def cut(func):
    def wrapper(**kwargs):
        if 'cut' in kwargs:
            pass
        else:
            return func(kwargs)
    return wrapper

def slider(func):
    def wrapper(**kwargs):
        fig = go.Figure()
        if 'slider' in kwargs:
            slider = kwargs['slider']
            key, args = list(slider['args'].items())[0]
            N_slider = len(args)
            for i in range(len(args)):
                traces = func(**{key: args[i]}, **kwargs)
                N_traces = 1
                if type(traces) == list:
                    N_traces = len(traces)
                    for trace in traces:
                        fig.add_trace(trace)
                else:
                    fig.add_trace(traces)
            steps = []
            for i in range(N_slider):
                step = dict(
                    method="update",
                    args=[{"visible": [False] * len(fig.data)}],  # layout attribute
                    label=str(i)
                )
                for j in range(N_traces):
                    step['args'][0]['visible'][N_traces*i+j] = True
                steps.append(step)
            if 'active' in slider:
                active = slider['active']
            else:
                active = 0
            for i in range(N_traces):
                fig.data[active+i].visible = True
            fig.update_layout(
                sliders = [dict(steps=steps, active=active)]
            )
            return fig
        else:
            return fig.add_trace(func(**kwargs))
    return wrapper