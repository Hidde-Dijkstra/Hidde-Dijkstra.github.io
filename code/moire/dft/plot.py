import numpy as np
from .decorators import change_directory, check
import plotly.graph_objects as go
from moire import lattice as dg
from moire.plot import figure, BandStructure
from scipy.interpolate import RectBivariateSpline

@change_directory('data_directory')
def plot_w90_bands(self, name='W90 bands', bs=None, α=1, **kwargs):
    if bs == None:
        bs = self.gen_bs()
        kwargs['reset_color'] = True
    unsorted_bands = np.load(self.prefix+'/w90_bands.npy')
    k_arr = np.load(self.prefix+'/w90_band_ticks.npy')
    bands = unsorted_bands[untangle(unsorted_bands, α=α)]
    bs.add_bands(bands, name=name, k_arr=k_arr, **kwargs)
    return bs

@check('qe')
@change_directory('data_directory')
def plot_qe_bands(self, name='QE bands', bs=None, α=1, **kwargs):
    if bs == None:
        bs = self.gen_bs()
        kwargs['reset_color'] = True
    unsorted_bands = np.load(self.prefix+'/qe_bands.npy')
    bands = unsorted_bands[untangle(unsorted_bands, α=α)]
    bs.add_bands(bands.T, name=name, **kwargs)
    return bs

@check('qe')
@change_directory('data_directory')
def plot_bands(self, α=1, dis_min=None, dis_max=None, print_N_wan=True, 
               color='blue', name='QE bands', bs=None, **kwargs):
    unsorted_bands = np.load(self.prefix+'/qe_bands.npy')
    bands = unsorted_bands[untangle(unsorted_bands, α=α)]
    if dis_max != None and dis_min != None:
        window = np.any(np.logical_and(dis_min<bands, bands<dis_max), axis=0)
        if print_N_wan:
            print('N wannier:', sum(window))
        self.band_plot.add_bands(bands[:, window], 'red', 'wannier bands', autocolor=False)
        self.band_plot.add_bands(bands[:, np.logical_not(window)], color, name, autocolor=False)
        self.band_plot.fig.add_hline(y=dis_min, line_dash="dash", line_width=1)
        self.band_plot.fig.add_hline(y=dis_max, line_dash="dash", line_width=1)
    else:
        self.band_plot.add_bands(bands.T, name=name, **kwargs)
    return self.band_plot.fig

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
def plot_wfc(self, dis_min=-np.inf, dis_max=np.inf, α=1, k_list=[], k_tags=[], 
             N_k=[], print_N_wan=True, ε_margin=1, **kwargs):
    band_dic = dict(
        hoverinfo='skip',
        line=dict(color='black', width=0.3))
    bs = BandStructure(cut=kwargs.get('cut'))
    bs.set_k_path(k_list, k_tags, N_k)
    unsorted_bands = np.load(self.prefix+'/projwfc/bands.npy')
    bands = unsorted_bands[untangle(unsorted_bands, α=α, spacing=bs.spacing)]
    window = np.any(np.logical_and(dis_min<bands, bands<dis_max), axis=0)
    if print_N_wan:
        print('N wannier:', sum(window))
    ε_range = [np.min(bands[:, window])-ε_margin, np.max(bands[:, window])+ε_margin]
    bs.add_bands(bands=bands[:, window].T, name='bands', autocolor=False, band_dic=band_dic)
    band_dic['line']['dash'] = 'dot'
    bs.add_bands(bands=bands[:, np.logical_not(window)].T, name='bands', 
                 autocolor=False, showlegend=False, plot_range=ε_range, 
                 band_dic=band_dic)
    return _plot_wfc(self=self, fig=bs.fig, α=α, dis_min=dis_min, dis_max=dis_max, 
                     bs=bs, **kwargs)

@change_directory('data_directory')
@figure
def _plot_wfc(
        self, dis_min=-np.inf, dis_max=np.inf, orbitals={}, 
        elements={}, α=1, size=15, spin=False, bs=BandStructure(), **kwargs):
    bands = np.load(self.prefix+'/projwfc/bands.npy')
    states = sort_states(np.load(self.prefix+'/projwfc/states.npy'), spin, orbitals, elements)
    occupations = np.load(self.prefix+'/projwfc/occupations.npy')# + 10**-6
    #occupations = occupations / np.expand_dims(np.sum(occupations, axis=2), axis=2)
    orbital_density = np.sum(occupations[:, :, states], axis=2)
    traces = []
    order = untangle(bands, α=α, spacing=bs.spacing)
    bands = bands[order]
    window = np.any(np.logical_and(dis_min<bands, bands<dis_max), axis=0)
    orbital_density = (orbital_density[order])[:, window]
    for i, band in enumerate(bands[:, window].T):
        marker_dic = dict(
            size=size*orbital_density[:, i],
            line=dict(width=1),
        )
        traces.append(go.Scattergl(
            x = bs.k_spacing, 
            y = band,
            opacity=0.7,
            text = orbital_density[:, i],
            mode = 'markers',
            marker = marker_dic
        ))
    return traces

@change_directory('data_directory')
def plot_wannier(self, orbital_dir=''):
    if orbital_dir != '':
        orbital_dir = '/' + orbital_dir
    origin = np.load(self.prefix+orbital_dir+'/wan_origins.npy')
    orbitals = np.load(self.prefix+orbital_dir+'/w90_orbitals.npy')
    vec_span = np.load(self.prefix+orbital_dir+'/w90_vec_span.npy')
    N_wan, n_1, n_2, n_3 = orbitals.shape
    n_x, n_y, n_z = convert_grid(orbitals[0, :, :, :]).shape
    n_ref = int(n_3/2) - int(n_2/2)
    X, Y, Z = np.mgrid[0:1:1j*n_x, 0:1:1j*n_y, n_ref/n_3:(n_ref+n_2)/n_3:1j*n_z]
    grid_origin = origin + [min([0, self.a[0][i], self.a[1][i]]) for i in range(3)]
    factors = sum([np.abs(vec_span[i]) for i in range(3)])
    iso_max = 0.9*np.min([np.amax(orbitals[i, :, :, :]) for i in range(N_wan)])
    x, y, z = [axis*factors[i]+grid_origin[i] for i, axis in enumerate([X, Y, Z])]
    W = np.max(factors[:2]/self.lattice.a)/2
    fig = self.lattice.plot(W=W, H=W, plot_3d=True)
    fig.update_layout(legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1),
        scene_aspectmode='cube'         
    )
    orbitals = np.roll(orbitals, int(n_3/2), axis=3)

    for i in range(self.w90_dic['num_wann']):
        fig.add_trace(go.Isosurface(
            x=x.flatten() - 2*sum([a[0] for a in self.a]),
            y=y.flatten(), #- sum([a[1] for a in self.a]),
            z=z.flatten() - self.Δz/2,
            value=convert_grid(orbitals[i, :, :, :]).flatten(),
            isomin=-iso_max,
            isomax=iso_max,
            opacity=0.6,
            surface_count=6,
            caps=dict(x_show=False, y_show=False, z_show=False),
            showlegend=True,
            name='orbital '+str(i+1)
        ))
    return fig   
    
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

def sort_states(states, spin, orbitals, elements):
    key_list, converted_states = convert_states(states, spin)
    wfc = {}
    for key in key_list:
        if key in orbitals:
            wfc[key] = list_transform(orbitals[key])
        elif key in elements:
            wfc[key] = list_transform(elements[key])
        else:
            wfc[key] = []
    return check_states(converted_states, wfc)

def convert_states(states, spin):
    if spin:
        atom, l, j, m_j = states.T
        m_l = m_j - np.sign((j-l)*m_j) / 2
        m_s = np.sign(j-l)*np.sign(m_j) / 2
        return ['atoms', 'l', 'm_l', 'm_s'], np.array([atom, l, m_l, m_s])
    else:
        atom, l, m = states.T
        return ['atoms', 'l', 'm'], np.array([atom, l, m])

def check_states(states, wfc):
    valid_states = []
    for i, state in enumerate(states.T):
        valid_state = True
        for key, value in zip(wfc.keys(), state):
            valid_state &= (value in wfc[key]) or (wfc[key] == [])
        if valid_state:
            valid_states.append(i)
    return valid_states

def list_transform(variable):
    if type(variable) == list:
        return variable
    else:
        return [variable]

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

def untangle(unsorted_bands, α=1, spacing=[]):
    M, N = unsorted_bands.shape
    if spacing == []:
        spacing = [M-1]
    points = [sum(spacing[:i]) for i in range(len(spacing)+1)]
    points[-1]
    bands = np.zeros((M, N))
    bands[:2] = unsorted_bands[:2]
    δ_bands = bands[1] - bands[0]
    index = np.zeros((M, N), dtype=int)
    for i, point in enumerate(points[:-1]):
        index[point] = range(N)
        index[point+1] = range(N)
        bands[point:point+2] = unsorted_bands[point:point+2]
        δ_bands = bands[point+1] - bands[point]
        for m in range(2+point, points[i+1]+1):
            δ_bands = α*(bands[m-1]-bands[m-2]) + (1-α)*δ_bands     
            bands_estimate = δ_bands + bands[m-1]
            bands = bands[:, np.argsort(bands_estimate)]
            δ_bands = δ_bands[np.argsort(bands_estimate)]
            index = index[:, np.argsort(bands_estimate)]
            bands[m] = unsorted_bands[m]
            index[m] = range(N)
    return np.arange(M)[:, None], index