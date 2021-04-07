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
def plot_wfc(self, i_ε):
    pass

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
