import os
from .decorators import check, time, change_directory
from moire.plot import BandStructure
import numpy as np

class DFT:

    from .write_file import write_qe, write_w90
    from .extract import extract_nscf, extract_w90, extract_projwfc, extract_orbitals
    from .lattice import import_lattice
    from .plot import plot_nscf, plot_wfc, plot_qe_bands, plot_w90_bands, plot_wannier

    a_1 = np.array([1, 0, 0])
    a_2 = np.array([-1/2, np.sqrt(3)/2, 0])
    a_3 = np.array([0, 0, 1])

    def __init__(self, prefix, work_dir, data_dir='', Δz=32):
        self.prefix = prefix
        self.work_directory = work_dir
        self.data_directory = data_dir
        self.qe = False
        self.w90 = False
        self.Δz = Δz

    def prepare_qe(self, qe_dic):
        self.qe = True
        self.qe_dic = {
            '&CONTROL': {
                'prefix': self.prefix,
                'outdir': './out',
                'verbosity': 'high'
            },
            '&SYSTEM': {
                'assume_isolated': '2D',
                'ibrav': 0,
                'nat': 0,
                'ntyp': 0
            },
            '&ELECTRONS': {},
            '&IONS': {}
        }
        for key in qe_dic.keys():
            for sub_key in qe_dic[key].keys():
                self.qe_dic[key][sub_key] = qe_dic[key][sub_key]
    
    @check('qe')
    def prepare_w90(self, w90_dic):
        self.w90 = True
        self.w90_dic = dict(guiding_centres=True, write_hr=True, wannier_plot=False,
                            bands_plot=True)
        if 'nbnd' in self.qe_dic['&SYSTEM']:
            self.w90_dic['num_bands'] = self.qe_dic['&SYSTEM']['nbnd']
        for key in w90_dic.keys():
            self.w90_dic[key] = w90_dic[key]
        self.orbital_centers = []
        a = self.qe_dic['&SYSTEM']['a']
        c = self.qe_dic['&SYSTEM']['c']
        self.a = [a*self.a_1, a*self.a_2, c*self.a_3]
        if self.w90_dic['wannier_plot']:
            self.w90_dic['wannier_plot_supercell']  = '3, 3, 1'

    @check('w90')
    def set_w90_window(self, win_min, win_max, froz_min, froz_max, num_iter):
        self.w90_dic['dis_win_min'] = win_min
        self.w90_dic['dis_win_max'] = win_max
        self.w90_dic['dis_froz_min'] = froz_min
        self.w90_dic['dis_froz_max'] = froz_max
        self.w90_dic['dis_num_iter'] = num_iter

    def run_qe(self, task, n_cores=4, nk=1):
        self._run('pw.x -nk '+str(nk), task, n_cores)

    @change_directory('work_directory')
    @time
    def run_w90(self, pp=True, n_cores=4):
        if pp:
            print('pp', os.system('wannier90.x -pp '+self.prefix+'.win'))
            self._run('pw2wannier90.x', 'pw2wan', n_cores)
        print('w90', os.system('wannier90.x '+self.prefix+'.win'))


    @change_directory('work_directory')
    @time
    def _run(self, executor, task, n_cores):
        command = (
            'mpirun -n ' + str(n_cores) + ' '
            + executor
            + ' -in ' 
            + self.prefix + '.' + task + '.in >'
            + self.prefix + '.' + task + '.out'
        )
        print(task, os.system(command))

    @change_directory('data_directory')
    def create_H(self, R_max, spin=False, print_N=True):
        orbitals = []
        for atom in self.lattice.unit_cell:
            for projection in atom.projections:
                orbitals.append(dict(position=atom.position, orbital=projection))
                if spin:
                    orbitals.append(dict(position=atom.position, orbital=projection))
        hop_list = np.load(self.prefix +'/w90_hopping.npy')
        hoppings = []
        N = 0
        for hop in hop_list:
            ΔL = orbitals[int(hop[3])-1]['position'] - orbitals[int(hop[4])-1]['position']
            ΔR = (hop[0]*self.a[0]+hop[1]*self.a[1]-ΔL*self.lattice.a)

            if np.dot(ΔR, ΔR) <= R_max**2:
                N += 1
                hop_mat = np.zeros((self.w90_dic['num_wann'], self.w90_dic['num_wann']), dtype=complex)
                hop_mat[int(hop[3])-1, int(hop[4])-1] = hop[5] + 1j*hop[6]
                hoppings.append([hop_mat, hop[0:2]])
        if print_N:
            print(N, 'hoppings found')
        return lambda k: sum([hop[0]*np.exp(2*np.pi*1j*(hop[1][0]*k[0]+hop[1][1]*k[1]))for hop in hoppings])

    def set_k_path(self, k_list, k_tags, N_k):
        self.k_list = k_list
        self.k_tags = k_tags
        self.N_k = N_k
        for i in range(len(k_list)):
            self.k_list[i] = np.array(k_list[i])
        bs = self.gen_bs()
        self.spacing = bs.spacing
        self.k_spacing = bs.k_spacing
    
    def set_k_grid(self, N):
        self.N_grid = N 
        self.k_grid = [np.array([i/N, j/N, 0]) for i in range(N) for j in range(N)] 

    def gen_bs(self):
        bs = BandStructure()
        bs.set_k_path(self.k_list, self.k_tags, self.N_k)
        return bs 