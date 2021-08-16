from .decorators import change_directory, save
from .io_processing import read_file, gen_lst
import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import griddata

class Wannier:

    def __init__(self, DFT, orbital_dir=''):
        self.DFT = DFT
        self.prefix = DFT.prefix
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.a_vec = DFT.W90.a_vec
        self.Δa = DFT.middle * sum(self.a_vec)/2
        self.dir = orbital_dir
        self.N = DFT.W90.w90_dic['num_wann']
    """
    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self):
        iso_list = []
        orbitals = {}
        for i in range(1, self.DFT.W90.w90_dic['num_wann']+1):
            f = read_file(f'{self.prefix}_{str(i).zfill(5)}.xsf')
            data = f.split('BEGIN_DATAGRID_3D_UNKNOWN\n')[1]
            data = data.split('\n')
            orbitals[f'{self.dir}/w90_N_grid'] = np.array([gen_lst(data[0], ' ', int)], 
                                                          dtype=int)
            orbitals[f'{self.dir}/wan_origins'] = gen_lst(data[1], ' ', float)
            orbitals[f'{self.dir}/w90_vec_span'] = [gen_lst(data[j], ' ', float) 
                                                    for j in range(2, 5)]
            iso_data = np.array([gen_lst(row, ' ', float) for row in data[5:-3]]).flatten()
            iso_data = iso_data.reshape(*np.flip(orbitals[f'{self.dir}/w90_N_grid']))
            iso_list.append(np.swapaxes(iso_data, 0, 2))
        orbitals[f'{self.dir}/w90_orbitals'] = np.array(iso_list)
        return orbitals
    """ 
  
    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self):
        origin, span, N, n_1, n_2, n_z = self.extract_iso_grid()
        xy_grid = [origin[:2]+np.dot(span, (i/n_1, j/n_2, 0))[:2]
                   for i in range(n_1) for j in range(n_2)]
        ΔR = np.max(span, axis=0)-np.min(span, axis=0) 
        x, y = [np.linspace(-min(ΔR)/2, min(ΔR/2), max(n_1, n_2)) + center
                for center in origin[:2]+np.dot(span, (1/2, 1/2, 0))[:2]]
        xi, yi = np.meshgrid(x, y)
        iso_vals = [[griddata(xy_grid, iso, (xi, yi), fill_value=0)
                     for iso in self.extract_iso(i, N, n_z, n_1*n_2)]
                    for i in range(1, self.N+1)]
        return {
            f'{self.dir}/iso_vals': iso_vals,
            f'{self.dir}/iso_x': x,
            f'{self.dir}/iso_y': y,
            f'{self.dir}/iso_z': [origin[2]+i/n_z*span[2][2] for i in range(n_z)]
        }
        
    def extract_iso_grid(self):
         with open(f'{self.prefix}_00001.xsf') as file:
            N = 0
            line = ''
            while line != 'BEGIN_DATAGRID_3D_UNKNOWN\n':
                line = next(file)
                N += 1
            n_1, n_2, n_z = gen_lst(next(file), ' ', int)
            origin =  gen_lst(next(file), ' ', float)
            span = [np.array(gen_lst(next(file), ' ', float)) 
                        for _ in range(3)]
            return origin, span, N, n_1, n_2, n_z
                
    def extract_iso(self, i, N, n_z, n_xy):
        with open(f'{self.prefix}_{str(i).zfill(5)}.xsf') as file:
            return np.array(
                [iso_val for row in file.read().split('\n')[N+5:-3]
                 for iso_val in gen_lst(row, ' ', float)]
                ).reshape((n_z, n_xy))
            
    @change_directory('data_directory')      
    def plot(self, W=3, step=1, visibility=True, N_iso=3, iso_factor=0.9) :
        x, y, z = [np.load(f'{self.prefix}/{self.dir}/iso_{axis}.npy')
                   for axis in ['x', 'y', 'z']]
        iso_vals = np.nan_to_num(np.load(f'{self.prefix}/{self.dir}/iso_vals.npy'))
        fig = self.DFT.lattice.plot(W=W, H=W, plot_3d=True)
        fig.update_layout(legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1),
            scene_aspectmode='cube'         
        )
        xcrop = np.where((np.mean(x)-W<x) * (x<np.mean(x)+W))[0]
        ycrop = np.where((np.mean(y)-W<y) * (y<np.mean(y)+W))[0]
        zcrop = np.where((np.mean(z)-W<z) * (z<np.mean(z)+W))[0]
        xi, yi, zi = np.meshgrid(x[xcrop[0]:xcrop[-1]], 
                                 y[ycrop[0]:ycrop[-1]], 
                                 z[zcrop[0]:zcrop[-1]])
        cropped_iso = np.swapaxes(iso_vals[::step,
                                           zcrop[0]:zcrop[-1],
                                           ycrop[0]:ycrop[-1], 
                                           xcrop[0]:xcrop[-1]],
                                  1, 3)
        iso_max = np.max(np.abs(cropped_iso))
        for i, iso_val in enumerate(cropped_iso):
            if i == 0:
                visible = True
            else:
                visible = visibility
            #"""
            fig.add_trace(go.Isosurface(
                x=xi.flatten() - np.mean(x),
                y=yi.flatten() - np.mean(y),
                z=zi.flatten() - np.mean(z),
                value=iso_val.flatten(),
                isomin=-iso_max*iso_factor,
                isomax=iso_max*iso_factor,
                opacity=0.6,
                surface_count=2*N_iso,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showlegend=True,
                visible=visible,
                name='orbital '+str(i+1)
            ))
            #"""
        return fig  
    
    """
    @change_directory('data_directory')
    def plot(self, step=1, visibility=True, N_iso=3, iso_factor=0.9):
        origin = np.load(f'{self.prefix}/{self.dir}/wan_origins.npy')
        orbitals = np.load(f'{self.prefix}/{self.dir}/w90_orbitals.npy')
        vec_span = np.load(f'{self.prefix}/{self.dir}/w90_vec_span.npy')
        N_wan, n_1, n_2, n_3 = orbitals.shape
        n_x, n_y, n_z = convert_grid(orbitals[0, :, :, :]).shape
        n_ref = int(n_3/2) - int(n_2/2)
        X, Y, Z = np.mgrid[0:1:1j*n_x, 0:1:1j*n_y, n_ref/n_3:(n_ref+n_2)/n_3:1j*n_z]
        grid_origin = origin + [min([0, self.a_vec[0][i], self.a_vec[1][i]]) 
                                for i in range(3)]
        factors = sum([np.abs(vec_span[i]) for i in range(3)])
        iso_max = np.min([np.amax(orbitals[i, :, :, :]) for i in range(N_wan)])
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
                isomin=-iso_max*iso_factor,
                isomax=iso_max*iso_factor,
                opacity=0.6,
                surface_count=2*N_iso,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showlegend=True,
                visible=visible,
                name='orbital '+str(i+1)
            ))
        return fig  
    """

def convert_grid(A):
    n_1, n_2, n_3 = A.shape
    grid = np.zeros((n_1+int(n_2/2)-1, int(n_2/2), n_3))
    for i in range(n_1):
        for j in range(int(n_2/2)):
            grid[int(n_2/2)-1+i-j, j, :] = A[i, 2*j, :]
    n = int(n_3/2)
    Δn = int(n_2/2)
    return grid[::2, :, n-Δn:n-Δn+n_2]