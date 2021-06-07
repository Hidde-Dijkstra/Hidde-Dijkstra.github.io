from .decorators import change_directory, save
from .io_processing import read_file, gen_lst
import numpy as np
import plotly.graph_objects as go

class Wannier:

    def __init__(self, DFT, orbital_dir=''):
        self.DFT = DFT
        self.prefix = DFT.prefix
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.a_vec = DFT.W90.a_vec
        if orbital_dir != '':
            self.dir = '/' + orbital_dir + '/'
        else:
            self.dir = '/'

    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self):
        iso_list = []
        orbital_data = {}
        for i in range(1, self.DFT.W90.w90_dic['num_wann']+1):
            f = read_file(self.prefix+'_'+str(i).zfill(5)+'.xsf')
            data = f.split('BEGIN_DATAGRID_3D_UNKNOWN\n')[1]
            data = data.split('\n')
            orbital_data[self.dir+'w90_N_grid'] = np.array([gen_lst(data[0], ' ', int)], dtype=int)
            orbital_data[self.dir+'wan_origins'] = gen_lst(data[1], ' ', float)
            orbital_data[self.dir+'w90_vec_span'] = [gen_lst(data[j], ' ', float) for j in range(2, 5)]
            iso_data = np.array([gen_lst(row, ' ', float) for row in data[5:-3]]).flatten()
            iso_data = iso_data.reshape(*np.flip(orbital_data[self.dir+'w90_N_grid']))
            iso_list.append(np.swapaxes(iso_data, 0, 2))
        orbital_data[self.dir+'w90_orbitals'] = np.array(iso_list)
        return orbital_data

    @change_directory('data_directory')
    def plot(self, step=1, visibility=True):
        origin = np.load(self.prefix+self.dir+'wan_origins.npy')
        orbitals = np.load(self.prefix+self.dir+'w90_orbitals.npy')
        vec_span = np.load(self.prefix+self.dir+'w90_vec_span.npy')
        N_wan, n_1, n_2, n_3 = orbitals.shape
        n_x, n_y, n_z = convert_grid(orbitals[0, :, :, :]).shape
        n_ref = int(n_3/2) - int(n_2/2)
        X, Y, Z = np.mgrid[0:1:1j*n_x, 0:1:1j*n_y, n_ref/n_3:(n_ref+n_2)/n_3:1j*n_z]
        grid_origin = origin + [min([0, self.a_vec[0][i], self.a_vec[1][i]]) for i in range(3)]
        factors = sum([np.abs(vec_span[i]) for i in range(3)])
        iso_max = 0.9*np.min([np.amax(orbitals[i, :, :, :]) for i in range(N_wan)])
        x, y, z = [axis*factors[i]+grid_origin[i] for i, axis in enumerate([X, Y, Z])]
        W = np.max(factors[:2]/self.DFT.lattice.a)/2
        fig = self.DFT.lattice.plot(W=W, H=W, plot_3d=True)
        fig.update_layout(legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1),
            scene_aspectmode='cube'         
        )
        orbitals = np.roll(orbitals, int(n_3/2), axis=3)

        for i in range(0, self.DFT.W90.w90_dic['num_wann'], step):
            if i == 0:
                visible = True
            else:
                visible = visibility
            fig.add_trace(go.Isosurface(
                x=x.flatten() - 2*sum([a[0] for a in self.a_vec]),
                y=y.flatten(), #- sum([a[1] for a in self.a]),
                z=z.flatten() - self.DFT.Δz/2,
                value=convert_grid(orbitals[i, :, :, :]).flatten(),
                isomin=-iso_max,
                isomax=iso_max,
                opacity=0.6,
                surface_count=6,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showlegend=True,
                visible=visible,
                name='orbital '+str(i+1)
            ))
        return fig  


def convert_grid(A):
    n_1, n_2, n_3 = A.shape
    grid = np.zeros((n_1+int(n_2/2)-1, int(n_2/2), n_3))
    for i in range(n_1):
        for j in range(int(n_2/2)):
            grid[int(n_2/2)-1+i-j, j, :] = A[i, 2*j, :]
    n = int(n_3/2)
    Δn = int(n_2/2)
    return grid[::2, :, n-Δn:n-Δn+n_2]