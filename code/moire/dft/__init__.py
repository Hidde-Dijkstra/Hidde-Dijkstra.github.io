import os
from .decorators import check, time, change_directory
from moire.bandstructure import BandStructure
import numpy as np

class DFT:

    from .write_file import write_qe, write_w90
    from .extract import extract_nscf, extract_w90_bands
    from .lattice import import_lattice
    from .plot import plot_nscf

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
        self.w90_dic = {
            'guiding_centres': 'true',
            'write_hr': 'true',
            'wannier_plot': 'true',
            'bands_plot': 'true',
            'wannier_plot_supercell': [3, 3, 1],
            'spinors': 'false'
        }
        if 'nbnd' in self.qe_dic['&SYSTEM']:
            self.w90_dic['num_bands'] = self.qe_dic['&SYSTEM']['nbnd']
        for key in w90_dic.keys():
            self.w90_dic[key] = w90_dic[key]
        self.orbital_centers = []

    def run_qe(self, task, n_cores=4, nk=1):
        self._run('pw.x -nk '+str(nk), task, n_cores)

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

    def set_k_path(self, k_list, k_tags, N_k):
        self.k_list = k_list
        self.k_tags = k_tags
        self.N_k = N_k
        for i in range(len(k_list)):
            k_list[i] = np.array(k_list[i])
        self.band_plot = BandStructure()
        self.band_plot.set_k_path(k_list, k_tags, N_k)
    
    def set_k_grid(self, N):
        self.N_grid = N 
        self.k_grid = [np.array([i/N, j/N, 0]) for i in range(N) for j in range(N)]  