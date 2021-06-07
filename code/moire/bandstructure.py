from .plot.figure import figure
from plotly.subplots import make_subplots
import numpy as np
import plotly.graph_objects as go


class BandStructure:
    
    def __init__(self, cut=None):
        self.fig = make_subplots(
            rows = 1 + (cut!=None), 
            cols = 1, 
            vertical_spacing = 0.05, 
            shared_xaxes=True,
            x_title='Momentum', 
            y_title= 'Energy (eV)',
        )
        self.fig.update_layout(legend= {'itemsizing': 'constant'})
        self.cut = cut
        self.spacing = None
    
    def plot_from_H(self, name='', showlegend=True, autocolor=True, **kwargs):
        if 'legend' not in kwargs:
            kwargs['legend'] = dict(args=dict(_=[None]), tags=[name], autocolor=autocolor)
        elif 'autocolor' not in kwargs['legend']:
            kwargs['legend']['autocolor'] = True
        self._plot_from_H(fig=self.fig, cut=self.cut, showlegend=showlegend, **kwargs)

    def add_bands(self, bands, name='', showlegend=True, autocolor=True, k_arr=None, **kwargs):
        if 'legend' not in kwargs:
            kwargs['legend'] = dict(args=dict(_=[None]), tags=[name], autocolor=autocolor)
        self._plot_bands(bands=bands, fig=self.fig, cut=self.cut, showlegend=showlegend,
            k_arr=k_arr, **kwargs)
    @figure
    def _plot_from_H(self, H=lambda k: 0, N_input=2, untangle_bands=True, α=1, **kwargs):
        N = len(H(np.zeros(N_input)))
        bands = np.zeros((N, len(self.k_path)))
        if 'band_dic' in kwargs:
            band_dic = kwargs['band_dic']
        else:
            band_dic = {}
        for j in range(len(self.k_path)):
            bands[:, j] = np.linalg.eigvalsh(H(self.k_path[j]))
        traces = []
        if untangle_bands:
            bands = bands[untangle(bands, α=α, spacing=self.spacing)]
        for band in bands:
            traces.append(go.Scatter(x=self.k_spacing, y=band, **band_dic))
        return traces

    @figure        
    def _plot_bands(self, bands=[], k_arr=None, untangle_bands=True, α=1, **kwargs):
        if 'band_dic' in kwargs:
            band_dic = kwargs['band_dic']
        else:
            band_dic = {}
        traces = []
        if np.all(k_arr == None):
            x = self.k_spacing
            spacing = self.spacing
        else:
            x = k_arr
            spacing = kwargs.get('spacing')
        if untangle_bands:
            bands = bands[untangle(bands, α=α, spacing=spacing)]
        for band in bands:
            traces.append(go.Scatter(x=x, y=band, **band_dic))
        return traces
    
    def set_k_path(self, k_list, k_tags, N_k):
        k_list = [np.array(k) for k in k_list]
        k_norms = [np.linalg.norm(k_list[i+1]-k_list[i]) for i in range(len(k_list)-1)]
        if type(N_k) == int:
            n = [int(N_k*k_norms[i]/sum(k_norms)) for i in range(len(k_norms))] 
            ΔN_k = N_k - sum(n)
            if ΔN_k != 0:
                rest = np.array(k_norms) - sum(k_norms) * np.array(n)/N_k
                np.array(n, dtype=int)[rest.argsort()[-ΔN_k:]] += 1
            self.spacing = n
        else:
            self.spacing  = N_k
        k_path = []
        k_spacing = []
        for i, Δn in enumerate(self.spacing):
            k_path += [k_list[i] + (k_list[i+1]-k_list[i])*j/Δn for j in range(Δn)]
            k_spacing += [sum(k_norms[:i]) + (k_norms[i])*j/Δn for j in range(Δn)]
        k_path.append(k_list[-1])
        k_spacing.append(sum(k_norms))
        self.k_spacing = np.array(k_spacing) / sum(k_norms)
        self.k_path = k_path
        tickvals = [self.k_spacing[sum(self.spacing[:i])] for i in range(len(self.spacing)+1)]
        self.fig.update_xaxes(ticktext=k_tags, tickvals=tickvals)


def untangle(unsorted_bands, α=1, spacing=None):
    bands = np.sort(unsorted_bands, axis=0)
    N, M = bands.shape
    if spacing == None:
        spacing = [M-1]
    points = [sum(spacing[:i]) for i in range(len(spacing)+1)]
    band_slices = np.zeros((M, N))
    index = np.zeros((M, N), dtype=int)
    for i, point in enumerate(points[:-1]):
        index[point] = range(N)
        index[point+1] = range(N)
        band_slices[point:point+2] = bands[:, point:point+2].T
        δ_bands = band_slices[point+1] - band_slices[point]
        for m in range(2+point, points[i+1]+1):
            δ_bands = α*(band_slices[m-1]-band_slices[m-2]) + (1-α)*δ_bands   
            bands_estimate = δ_bands + band_slices[m-1]
            band_slices = band_slices[:, np.argsort(bands_estimate)]
            δ_bands = δ_bands[np.argsort(bands_estimate)]
            index = index[:, np.argsort(bands_estimate)]
            band_slices[m] = bands[:, m]
            index[m] = range(N)
    return index.T, np.arange(M)[:, None].T