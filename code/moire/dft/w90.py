from moire.plot.heatmap import gen_heatmap
from .decorators import change_directory, save
from .io_processing import read_file, gen_lst, write_file, join_grid, join_grid_point
from collections import defaultdict
import numpy as np
import os

class W90:

    def __init__(self, DFT, w90_dic):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.w90_dic = dict(guiding_centres=True, write_hr=True, wannier_plot=False,
                            bands_plot=True, spinors=False)
        if 'nbnd' in DFT.QE.qe_dic['&SYSTEM']:
            self.w90_dic['num_bands'] = DFT.QE.qe_dic['&SYSTEM']['nbnd']
        for key in w90_dic.keys():
            self.w90_dic[key] = w90_dic[key]
        self.orbital_centers = []
        self.a = DFT.QE.qe_dic['&SYSTEM']['a']
        c = DFT.QE.qe_dic['&SYSTEM']['c']
        self.a_vec = [self.a*DFT.a_1, self.a*DFT.a_2, self.a*DFT.a_3]
        if self.w90_dic['wannier_plot']:
            self.w90_dic['wannier_plot_supercell']  = '3, 3, 1'

    def set_window(self, win_min, win_max, froz_min, froz_max, num_iter):
        self.w90_dic['dis_win_min'] = win_min
        self.w90_dic['dis_win_max'] = win_max
        self.w90_dic['dis_froz_min'] = froz_min
        self.w90_dic['dis_froz_max'] = froz_max
        self.w90_dic['dis_num_iter'] = num_iter

    @change_directory('work_directory')
    def write(self, pw2wan=True):
        if pw2wan:
            self.DFT.QE.write_pw2wan(self.w90_dic['wannier_plot'])
        w90_file = ''
        for setting, value in self.w90_dic.items():
            if type(value) == list:
                w90_file += setting + ' = ' + ', '.join([str(item) for item in value]) + '\n'
            else:
                w90_file += setting + ' = ' + str(value).lower() + '\n'
        w90_file += '\nbegin unit_cell_cart\n'
        for a_i in self.a_vec:
            w90_file += '  ' + join_grid_point(a_i) + '\n'
        w90_file += 'end unit_cell_cart\n\nBegin projections'
        for atom_name, atom in self.DFT.atoms.items():
            projections = atom.get('projections')
            if projections != [] and projections != None:
                w90_file += '\n  ' + atom_name +':  ' + '; '.join(projections)
        w90_file += '\nEnd projections\n\nBegin atoms_cart\nang\n'
        for atom_name , atom in self.DFT.atoms.items():
            for site in atom['loc']:
                w90_file += '  ' + atom_name + '   ' + join_grid_point(site) + '\n'
        w90_file += 'End atoms_cart\n\nBegin Kpoint_Path\n'
        for i in range(len(self.DFT.k_tags)-1):
            w90_file += '  ' + self.DFT.k_tags[i] + ' ' + join_grid_point(self.DFT.k_list[i]) + '   '
            w90_file += self.DFT.k_tags[i+1] + ' ' + join_grid_point(self.DFT.k_list[i+1]) + '\n'
        w90_file += 'End Kpoint_Path\n\nmp_grid = ' + str(self.DFT.N_grid) +', ' + str(self.DFT.N_grid) +', 1'
        w90_file += '\n\nBegin kpoints\n' + join_grid(self.DFT.k_grid, weight=False) + 'End kpoints'
        write_file(self.DFT.prefix+'.win', w90_file)  

    @change_directory('work_directory')
    def run(self, pp=True, n_cores=4):
        if pp:
            print('pp', os.system('wannier90.x -pp '+self.prefix+'.win'))
            self.DFT.QE._run('pw2wannier90.x', 'pw2wan', n_cores)
        print('w90', os.system('wannier90.x '+self.prefix+'.win'))

    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self):
        band_file = read_file(self.prefix+'_band.dat')
        w90_data = [[gen_lst(elem, ' ', float) for elem in gen_lst(lst, '\n')] for lst in gen_lst(band_file, '\n  \n') ]
        w90_bands = np.array(w90_data)[:, :, 1]
        w90_k = np.array(w90_data)[0, :, 0]

        ticks_file = read_file(self.prefix+'_band.gnu')
        x_ticks = [elem.split("  ") for elem in (ticks_file.split("(")[1]).split(")")[0].split(",")]
        k_vals = [float(elem[1])-0.0001 for elem in x_ticks]
        n_ticks = [np.where(w90_k >= k)[0][0] for k in k_vals[:-1]]
        n_ticks += [len(w90_k)-1]
        x_arr = []
        ticks = [sum(self.DFT.spacing[:i]) for i in range(len(self.DFT.spacing)+1)]
        for i in range(len(ticks)-1):
            Δx = self.DFT.k_spacing[ticks[i+1]] - self.DFT.k_spacing[ticks[i]]
            Δn = n_ticks[i+1] - n_ticks[i]
            x_arr += [n*Δx/Δn+self.DFT.k_spacing[ticks[i]] for n in range(Δn)]
        x_arr += [self.DFT.k_spacing[-1]]
        hopping_file = read_file(self.DFT.prefix+'_hr.dat')
        hop_list = [gen_lst(lst, '   ', float) for lst in gen_lst(hopping_file, '\n')[8:-1]]    
        return dict(w90_bands=w90_bands, w90_band_ticks=x_arr, w90_hopping=hop_list)

    @change_directory('data_directory')
    def plot_bands(self, name='W90 bands', bs=None, **kwargs):
        if bs == None:
            bs = self.DFT.gen_bs()
            kwargs['reset_color'] = True
        bands = np.load(self.prefix+'/w90_bands.npy')
        k_arr = np.load(self.prefix+'/w90_band_ticks.npy')
        bs.add_bands(bands, name=name, k_arr=k_arr, **kwargs)
        return bs

    @change_directory('data_directory')
    def plot_potential(self, spin=False, atom_list=None):
        orbitals = self.get_orbitals()
        hops = np.load(self.prefix +'/w90_hopping.npy')       
        on_site = (hops[:, 0]==0) & (hops[:, 1]==0) & (hops[:, 3]==hops[:, 4])
        hop_dic = defaultdict(list)
        for hop in hops[on_site]:
            orbital = orbitals[int(hop[3])-1]
            key = orbital['orbital'] + (' '+str(orbital['spin'])) * spin
            if atom_list == None or orbital['i'] in atom_list:
                hop_dic[key].append(np.array([*orbital['position'][:2], np.abs(hop[5] + 1j*hop[6])]))
        slider = dict(args=dict(positions=list(hop_dic.values())), tags=list(hop_dic.keys()))
        return gen_heatmap(slider=slider, a=self.a)

    @change_directory('data_directory')
    def create_H(self, R_max, print_N=True):
        orbitals = self.get_orbitals()
        hop_list = np.load(self.prefix +'/w90_hopping.npy')
        hoppings = []
        N = 0
        for hop in hop_list:
            ΔL = orbitals[int(hop[3])-1]['position'] - orbitals[int(hop[4])-1]['position']
            ΔR = (hop[0]*self.a_vec[0]+hop[1]*self.a_vec[1]-ΔL*self.DFT.lattice.a)

            if np.dot(ΔR, ΔR) <= R_max**2:
                N += 1
                hop_mat = np.zeros((self.w90_dic['num_wann'], self.w90_dic['num_wann']), dtype=complex)
                hop_mat[int(hop[3])-1, int(hop[4])-1] = hop[5] + 1j*hop[6]
                hoppings.append([hop_mat, hop[0:2]])
        if print_N:
            print(N, 'hoppings found')
        return lambda k: sum([hop[0]*np.exp(2*np.pi*1j*(hop[1][0]*k[0]+hop[1][1]*k[1]))for hop in hoppings])

    def get_orbitals(self):
        orbitals = []
        for i, atom in enumerate(self.DFT.lattice.unit_cell):
            for projection in atom.projections:
                orbitals.append(dict(position=atom.position, orbital=projection, spin=-1, i=i))
                if self.w90_dic['spinors']:
                    orbitals.append(dict(position=atom.position, orbital=projection, spin=1, i=i))
        return orbitals