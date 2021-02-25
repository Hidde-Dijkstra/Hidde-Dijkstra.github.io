# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.8.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Numerics

# +
import numpy as np
import plotly.graph_objects as go
import os
from datetime import datetime

from moire import bandstructure as bs
from moire import diagram_gen as dg


# -

class DFT:
    
    qe = False
    w90 = False
    
    def __init__(self, directory, prefix, lattice_vectors, data_dir='dft_data', Δz=32):
        self.prefix = prefix
        self.dir = directory
        self.data_dir = data_dir
        self.dir_exists(data_dir)
        self.dir_exists(data_dir+'/'+prefix)
        self.cwd = os.getcwd()
        self.a = lattice_vectors
        self.Δz = Δz
        
    def dir_exists(self, data_dir):
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
            print('directory created:', data_dir)
            
    def set_atoms(self, atoms):
        self.atoms = atoms
        for atom in atoms.keys():
            for j in range(len(atoms[atom]['loc'])):
                self.atoms[atom]['loc'][j][2] += self.Δz/2
        if self.qe:
            self.qe_dic['&SYSTEM']['nat'] = sum([len(atom['loc']) for atom in atoms.values()])
            self.qe_dic['&SYSTEM']['ntyp'] = len(atoms)
        if self.w90:
            projections = '_'.join([key+'_'+'_'.join(val['projections']) for key, val in atoms.items() if 'projections' in val])
            self.orbital_dir = self.data_dir + '/' + self.prefix + '/' + projections
            for atom in atoms.values():
                if 'projections' in atom:
                    for i in range(len(atom['loc'])):
                        for _ in range(len(atom['projections'])*(1+(self.w90_dic['spinors']=='true'))):
                            self.orbital_centers.append(atom['loc'][i])
            if len(self.orbital_centers) != self.w90_dic['num_wann']:
                raise Exception('Number of wannier orbitals and number of projections do not match: '+
                                str(len(self.orbital_centers))+' != '+str(self.w90_dic['num_wann']))
                    
    def set_k_path(self, k_list, k_tags, N_k):
        self.k_list = k_list
        self.k_tags = k_tags
        self.N_k = N_k
        for i in range(len(k_list)):
            k_list[i] = np.array(k_list[i])
        self.band_plot = bs.BandStructure()
        self.band_plot.set_k_path(k_list, k_tags, N_k)
    
    def set_k_grid(self, N):
        self.N_grid = N 
        self.k_grid = [np.array([i/N, j/N, 0]) for i in range(N) for j in range(N)]
        
    def join_grid_point(self, grid_point):
        return '   '.join(['{:.9f}'.format(item)[0:9] for item in grid_point])
        
    def join_k_grid(self, k_grid, weight=True):
        k_str = ''
        for grid_point in k_grid:
            k_str += '  ' + self.join_grid_point(grid_point) + ' 1'*weight+'\n'
        return k_str
    
    def plot_bands(self, ticks=None):
        qe_path = self.data_dir+'/'+self.prefix +'/qe_bands.npy'
        w90_path = self.data_dir+'/'+self.prefix +'/w90_bands.npy'
        w90_ticks = self.data_dir+'/'+self.prefix +'/w90_band_ticks.npy'
        if os.path.exists(qe_path) and self.qe:
            self.band_plot.add_bands(np.load(qe_path), 'blue', 'qe bands')
        if os.path.exists(w90_path) and self.w90:
            self.band_plot.fig.update_yaxes(range=[self.w90_dic['dis_win_min'], self.w90_dic['dis_win_max']])
            self.band_plot.add_bands(np.load(w90_path), 'green', 'w90 bands', k_arr=np.load(w90_ticks))


class QE(DFT):
    
    qe = True
    
    def __init__(self, qe_dic, directory, prefix, lattice_vectors, data_dir='dft_data', Δz=32):
        super().__init__(directory, prefix, lattice_vectors, data_dir=data_dir, Δz=Δz)
        self.qe_dic = {
            '&CONTROL': {
                'prefix': prefix,
                'outdir': './out',
                'verbosity': 'high'
            },
            '&SYSTEM': {
                'assume_isolated': '2D',
                'ibrav': 0,
                'nat': 0,
                'ntyp': 0
            },
            '&ELECTRONS': {}
        }
        for key in qe_dic.keys():
            for sub_key in qe_dic[key].keys():
                self.qe_dic[key][sub_key] = qe_dic[key][sub_key]
                
    def write_qe_file(self, qe_type):
        file = ''
        self.qe_dic['&CONTROL']['calculation'] = qe_type
        for elem in self.qe_dic.keys():
            file += elem + '\n'
            for sub_elem in self.qe_dic[elem].keys():
                var = self.qe_dic[elem][sub_elem]
                if type(var) == str:
                    var = "'" + var + "'"
                elif type(var) == bool:
                    var =str(var).lower()
                else:
                    var = str(var)
                file += '  ' + sub_elem + ' = ' + var + '\n'
            file += '/\n'
        file += 'CELL_PARAMETERS angstrom\n'
        for a_i in self.a + [(0, 0, self.Δz)]:
            file += '  ' + self.join_grid_point(a_i) + '\n'
        file += 'ATOMIC_SPECIES\n'
        for atom in self.atoms.keys():
            file += '  ' + '   '.join([atom, str(self.atoms[atom]['weight']), "'"+self.atoms[atom]['pseudo']+"'"]) + '\n'
        file += 'ATOMIC_POSITIONS angstrom\n'
        for atom in self.atoms.keys():
            for site in self.atoms[atom]['loc']:
                file += '  ' + atom + '   ' + self.join_grid_point(site) + '\n'
        if qe_type == 'scf':
            file += 'K_POINTS automatic\n  9 9 1 1 1 1'
        elif qe_type == 'nscf':
            file +=  'K_POINTS crystal\n' + str(len(self.k_grid)) + '\n'
            file += self.join_k_grid(self.k_grid)
        elif qe_type == 'bands':
            file +=  'K_POINTS crystal\n' + str(len(self.band_plot.k_path)) + '\n'
            file += self.join_k_grid(self.band_plot.k_path)
        f = open(self.dir+'/'+self.prefix+'.'+qe_type+'.in', "w")
        f.write(file)
        f.close()
        
    def run_qe(self, qe_type, cores=4):
        os.chdir(self.dir)
        file = self.dir+'/'+self.prefix+'.'+qe_type
        startTime = datetime.now()
        print(qe_type, os.system('mpirun -n '+str(cores)+' pw.x < '+file+'.in > '+file+'.out'))
        print(datetime.now()-startTime)
        os.chdir(self.cwd)
        if qe_type == 'bands':
            self.extract_qe_bands()
    
    def extract_qe_bands(self):
        f = open(self.dir+'/'+self.prefix+'.bands.out', "r")
        N_data, band_data = f.read().split("End of band structure calculation")
        f.close()
        N_b = int(N_data.split("number of Kohn-Sham states=")[1].split("\n")[0])
        band_data = band_data.split("k = ")
        qe_bands = np.zeros((self.N_k, N_b))
        for i in range(0, self.N_k):
            k, band_list = band_data[i+1].split("bands (ev):\n\n")
            if i == self.N_k - 1:
                band_list = band_list.split("Writing")[0]
            qe_bands[i, :] = [float(band) for band in band_list.replace("\n", "").split(" ") if band != ""]
        np.save(self.data_dir+'/'+self.prefix +'/qe_bands', qe_bands)


class W90(QE):
    
    w90 = True
    
    def __init__(self, w90_dic, qe_dic, directory, prefix, lattice_vectors, data_dir='dft_data', Δz=32):
        super().__init__(qe_dic, directory, prefix, lattice_vectors, data_dir=data_dir, Δz=Δz)
        self.w90_dic = {
            'guiding_centres': 'true',
            'write_hr': 'true',
            'wannier_plot': 'true',
            'bands_plot': 'true',
            'wannier_plot_supercell': [3, 3, 1],
            'spinors': 'false'
        }
        if 'nbnd' in qe_dic['&SYSTEM']:
            self.w90_dic['num_bands'] = qe_dic['&SYSTEM']['nbnd']
        for key in w90_dic.keys():
            self.w90_dic[key] = w90_dic[key]
        self.orbital_centers = []
            
    def set_window(self, win_min, win_max, froz_min, froz_max, num_iter):
        self.w90_dic['dis_win_min'] = win_min
        self.w90_dic['dis_win_max'] = win_max
        self.w90_dic['dis_froz_min'] = froz_min
        self.w90_dic['dis_froz_max'] = froz_max
        self.w90_dic['dis_num_iter'] = num_iter
        
    def write_w90_file(self):
        self.dir_exists(self.orbital_dir)
        os.chdir(self.dir)
        file = "&inputpp\n   outdir = '" + self.qe_dic['&CONTROL']['outdir'] +"'\n   "
        file += "prefix = '" + self.qe_dic['&CONTROL']['prefix'] +"'\n   "
        file += "seedname = '" + self.qe_dic['&CONTROL']['prefix'] +"'\n   "
        file += 'write_unk = .true.\n   write_mmn = .true.\n   write_amn = .true.\n/\n'
        f = open(self.dir+'/'+self.qe_dic['&CONTROL']['prefix']+'.pw2wan', 'w')
        f.write(file)
        f.close()
        file = ''
        for key in self.w90_dic.keys():
            var = self.w90_dic[key]
            if type(var) == list:
                file += key + ' = ' + ', '.join([str(item) for item in var]) + '\n'
            else:
                file += key + ' = ' + str(var) + '\n'
        file += '\nbegin unit_cell_cart\n'
        for a_i in self.a + [(0, 0, self.Δz)]:
            file += '  ' + self.join_grid_point(a_i) + '\n'
        file += 'end unit_cell_cart\n\nBegin projections'
        for atom in self.atoms.keys():
            if 'projections' in self.atoms[atom]:
                file += '\n  ' + atom +':  ' + '; '.join(self.atoms[atom]['projections'])
        file += '\nEnd projections\n\nBegin atoms_cart\nang\n'
        for atom in self.atoms.keys():
            for site in self.atoms[atom]['loc']:
                file += '  ' + atom + '   ' + self.join_grid_point(site) + '\n'
        file += 'End atoms_cart\n\nBegin Kpoint_Path\n'
        for i in range(len(self.k_tags)-1):
            file += '  ' + self.k_tags[i] + ' ' + self.join_grid_point(self.k_list[i]) + '   '
            file += self.k_tags[i+1] + ' ' + self.join_grid_point(self.k_list[i+1]) + '\n'
        file += 'End Kpoint_Path\n\nmp_grid = ' + str(self.N_grid) +', ' + str(self.N_grid) +', 1'
        file += '\n\nBegin kpoints\n' + self.join_k_grid(self.k_grid, weight=False) + 'End kpoints'
        f = open(self.dir+'/'+self.prefix+'.win', "w")
        f.write(file)
        f.close()
        os.chdir(self.cwd)
        
    def run_w90(self, pp=False, cores=4):
        os.chdir(self.dir)
        file = self.dir+'/'+self.prefix
        if pp:
            startTime = datetime.now()
            print('pp', os.system('wannier90.x -pp '+file))
            print(datetime.now()-startTime)
            startTime = datetime.now()
            print('pw2wan', os.system('mpirun -n '+str(cores)+' pw2wannier90.x < '+file+'.pw2wan > '+self.dir+'/pw2wan.out'))
            print(datetime.now()-startTime)
        startTime = datetime.now()
        print('w90', os.system('wannier90.x '+file))
        print(datetime.now()-startTime)
        os.chdir(self.cwd)
        f = open(self.dir+'/'+self.prefix+'_band.dat', "r")
        w90_bands = np.array([[[float(x) for x in elem.split(" ") if x!=""] for elem in lst.split("\n") if elem!=""] for lst in f.read().split("\n  \n") if lst!=""])
        f.close()
        f = open(self.dir+'/'+self.prefix+'_hr.dat', "r")
        hop_list = np.array([[float(item) for item in lst.split("   ") if item!=''] for lst in f.read().split("\n")[8:-1]])
        f.close()
        np.save(self.data_dir+'/'+self.prefix +'/hopping_elements', hop_list)
        self.extract_w90_bands()
        self.extract_orbitals()
        
    def extract_w90_bands(self):
        f = open(self.dir+'/'+self.prefix+'_band.dat', "r")
        w90_bands = np.array([[[float(x) for x in elem.split(" ") if x!=""] for elem in lst.split("\n") if elem!=""] for lst in f.read().split("\n  \n") if lst!=""])
        my_shape = np.shape(w90_bands)
        f.close()
        f = open(self.dir+'/'+self.prefix+'_band.gnu', "r")
        x_ticks = [elem.split("  ") for elem in (f.read().split("(")[1]).split(")")[0].split(",")]
        f.close()
        k_tags = [elem[0].split('"')[1] for elem in x_ticks]
        k_vals = [float(elem[1])-0.0001 for elem in x_ticks]
        n_ticks = [np.where(w90_bands[0, :, 0] >= k)[0][0] for k in k_vals[:-1]]
        n_ticks += [len(w90_bands[0, :, 0])-1]
        x_arr = []
        ticks = [sum(self.band_plot.spacing[:i]) for i in range(len(self.band_plot.spacing)+1)]
        for i in range(len(ticks)-1):
            Δx = ticks[i+1] - ticks[i]
            Δn = n_ticks[i+1] - n_ticks[i]
            x_arr += [n*Δx/Δn+ticks[i] for n in range(Δn)]
        x_arr += [self.N_k-1]
        np.save(self.data_dir+'/'+self.prefix +'/w90_bands', w90_bands[:, :, 1].T)
        np.save(self.data_dir+'/'+self.prefix +'/w90_band_ticks', np.array(x_arr))
        
    def extract_orbitals(self):
        iso_list = []
        for i in range(1, self.w90_dic['num_wann']+1):
            f = open(self.dir+'/'+self.prefix+'_'+str(i).zfill(5)+'.xsf', "r")
            data = f.read().split('BEGIN_DATAGRID_3D_UNKNOWN\n')[1]
            f.close()
            data = data.split('\n')
            N_grid = np.array([int(item) for item in data[0].split(' ') if item!=''], dtype=int)
            wan_origin = [float(item) for item in data[1].split(' ') if item!='']
            vec_span = np.array([[float(item) for item in data[i+2].split(' ') if item!=''] for i in range(3)])
            if i==1:
                np.save(self.data_dir+'/'+self.prefix +'/N_w90_grid_points', N_grid)
                np.save(self.data_dir+'/'+self.prefix +'/w90_grid_origin', wan_origin)
                np.save(self.data_dir+'/'+self.prefix +'/w90_vec_span', vec_span)
            iso_vals = np.array([float(item) for row in data[5:-3] for item in row.split(' ') if item!='']).reshape(*reversed(N_grid))
            iso_list.append(np.swapaxes(iso_vals, 0, 2))
        np.save(self.orbital_dir+'/w90_orbitals', np.array(iso_list))
        
    def create_H(self, R_max, show_amount_NN=False):
        hop_list = np.load(self.data_dir+'/'+self.prefix +'/hopping_elements.npy')
        vec_span = np.load(self.data_dir+'/'+self.prefix +'/w90_vec_span.npy')
        hoppings = []
        N = 0
        for hop in hop_list:
            ΔL = self.orbital_centers[int(hop[3])-1] - self.orbital_centers[int(hop[4])-1]
            ΔR = (hop[0]*self.a[0]+hop[1]*self.a[1]-ΔL)
            if np.dot(ΔR, ΔR) <= R_max**2:
                N += 1
                hop_mat = np.zeros((self.w90_dic['num_wann'], self.w90_dic['num_wann']), dtype=complex)
                hop_mat[int(hop[3])-1, int(hop[4])-1] = hop[5] + 1j*hop[6]
                hoppings.append([hop_mat, hop[0:2]])
        if show_amount_NN:
            print(N, 'hoppings found')
        return lambda k: sum([hop[0]*np.exp(2*np.pi*1j*(hop[1][0]*k[0]+hop[1][1]*k[1]))for hop in hoppings])
    
    def plot_tight_binding(self, R_max, name='tb', color='red', show_amount_NN=False):
        self.band_plot.plot_from_H(self.create_H(R_max, show_amount_NN=show_amount_NN), name, color, N_input=3)
            
    def convert_grid(self, A):
        n_1, n_2, n_3 = A.shape
        grid = np.zeros((n_1+int(n_2/2)-1, int(n_2/2), n_3))
        for i in range(n_1):
            for j in range(int(n_2/2)):
                grid[int(n_2/2)-1+i-j, j, :] = A[i, 2*j, :]
        n = int(n_3/2)
        Δn = int(n_2/2)
        return grid[::2, :, n-Δn:n-Δn+n_2]
    
    def plot_wannier(self, orbital_list, NN_distance, display_bool=True):
        origin = np.load(self.data_dir+'/'+self.prefix +'/w90_grid_origin.npy')
        orbitals = np.load(self.orbital_dir+'/w90_orbitals.npy')
        vec_span = np.load(self.data_dir+'/'+self.prefix +'/w90_vec_span.npy')
        N_wan, n_1, n_2, n_3 = orbitals.shape
        n_x, n_y, n_z = self.convert_grid(orbitals[0, :, :, :]).shape
        n_ref = int(n_3/2) - int(n_2/2)
        X, Y, Z = np.mgrid[0:1:1j*n_x, 0:1:1j*n_y, n_ref/n_3:(n_ref+n_2)/n_3:1j*n_z]
        grid_origin = origin + [min([0, self.a[0][i], self.a[1][i]]) for i in range(3)]
        factors = sum([np.abs(vec_span[i]) for i in range(3)])
        iso_max = 0.9*np.min([np.amax(orbitals[i, :, :, :]) for i in range(N_wan)])
        x, y, z = [axis*factors[i]+grid_origin[i] for i, axis in enumerate([X, Y, Z])]
        W = np.max(factors[:2])
        lattice = dg.Lattice3d(*self.a, W, W)
        fig = go.Figure()
        for key, atom in self.atoms.items():
            for loc in atom['loc']:
                lattice.add_atom(dg.LatticeAtom(loc, name=key, atom_color=atom['color']))
        lattice.gen_NN(NN_distance)
        lattice.draw_3d(fig, sum(self.a))
        fig.update_layout(legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1),
            scene_aspectmode='cube'         
        )
        if display_bool == True:
            display_bool = [True for i in orbital_list]
        for i, j_orb in enumerate(orbital_list):
            if display_bool[i]:
                visible=True
            else:
                visible='legendonly'
            fig.add_trace(go.Isosurface(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                value=self.convert_grid(orbitals[j_orb-1, :, :, :]).flatten(),
                isomin=-iso_max,
                isomax=iso_max,
                opacity=0.6,
                surface_count=6,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showlegend=True,
                visible=visible,
                name='orbital '+str(i)
            ))
        return fig
