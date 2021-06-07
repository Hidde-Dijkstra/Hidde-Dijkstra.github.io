from .decorators import change_directory, save, time
from .io_processing import read_file, gen_lst, write_file, join_grid, join_grid_point
import numpy as np
import os

class QE:

    def __init__(self, DFT, qe_dic, lattice):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.qe_dic = {
            '&CONTROL': {
                'prefix': DFT.prefix,
                'outdir': './out',
                'verbosity': 'high'
            },
            '&SYSTEM': {
                'assume_isolated': '2D',
                'ibrav': 4,
                'nat': len(lattice.unit_cell),
                'ntyp': len(lattice.atom_types),
                'a': lattice.a,
                'c': DFT.Δz
            },
            '&ELECTRONS': {},
            '&IONS': {}
        }
        for key in qe_dic.keys():
            for sub_key in qe_dic[key].keys():
                self.qe_dic[key][sub_key] = qe_dic[key][sub_key]

    @change_directory('work_directory')
    def write(self, qe_type, k_auto=9):
        qe_file = ''
        self.qe_dic['&CONTROL']['calculation'] = qe_type
        for section_name, section in self.qe_dic.items():
            qe_file += section_name + '\n'
            for setting, value in section.items():
                if type(value) == str:
                    value = "'" + value + "'"
                elif type(value) == bool:
                    value = str(value).lower()
                else:
                    value = str(value)
                qe_file += '  ' + setting + ' = ' + value + '\n'
            qe_file += '/\n'
        if self.qe_dic['&SYSTEM']['ibrav'] == 0:
            qe_file += 'CELL_PARAMETERS angstrom\n'
            for a_i in self.DFT.a + [(0, 0, self.DFT.Δz)]:
                qe_file += '  ' + join_grid_point(a_i) + '\n'
        qe_file += 'ATOMIC_SPECIES\n'
        for atom_name, atom in self.DFT.atoms.items():
            qe_file += '  ' + '   '.join([
                atom_name.ljust(2), 
                str(atom['weight']).ljust(6, '0'), 
                "'"+atom['pseudo_potential']+"'\n"
            ]) 
        qe_file += 'ATOMIC_POSITIONS angstrom\n'
        for atom_name, atom in self.DFT.atoms.items():
            for site in atom['loc']:
                qe_file += '  ' + atom_name.ljust(2) + '   ' + join_grid_point(site) + '\n'
        if qe_type == 'scf' or qe_type == 'relax':
            qe_file += 'K_POINTS automatic\n  '+str(k_auto)+' '+str(k_auto)+' 1 1 1 1'
        elif qe_type == 'nscf':
            qe_file +=  'K_POINTS crystal\n' + str(len(self.DFT.k_grid)) + '\n'
            qe_file += join_grid(self.DFT.k_grid)
        elif qe_type == 'bands':
            qe_file +=  'K_POINTS crystal\n' + str(len(self.DFT.k_path)) + '\n'
            qe_file += join_grid(self.DFT.k_path)
        write_file(self.prefix+'.'+qe_type+'.in', qe_file)


    def write_pw2wan(self, wannier_plot):
        control_dic = self.qe_dic['&CONTROL']
        pw2wan_file = "&inputpp\n   outdir = '" + control_dic['outdir'] +"'\n   "
        pw2wan_file += "prefix = '" + control_dic['prefix'] +"'\n   "
        pw2wan_file += "seedname = '" + control_dic['prefix'] +"'\n   "
        if wannier_plot:
            pw2wan_file += 'write_unk = .true.\n   '
        pw2wan_file  += 'write_mmn = .true.\n   write_amn = .true.\n/\n'
        write_file(control_dic['prefix']+'.pw2wan.in', pw2wan_file)

    def run(self, task, n_cores=4, nk=1):
        self._run('pw.x -nk '+str(nk), task, n_cores)
    
    @change_directory('work_directory')
    @time
    def _run(self, executor, task, n_cores):
        command = (
            'mpirun -n ' + str(n_cores) + ' '
            + executor
            + ' -in ' 
            + self.DFT.prefix + '.' + task + '.in >'
            + self.DFT.prefix + '.' + task + '.out'
        )
        print(task, os.system(command))

    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self, task='nscf'):
        nscf_file = read_file(self.prefix+'.'+task+'.out')
        nscf_data = nscf_file.split('End of band structure calculation')[1]
        nscf_data = gen_lst(nscf_data, 'k =')
        qe_nscf = []
        for line in nscf_data:
            nscf_point = line.split('bands (ev):\n\n')[1].split('occupation')[0]
            nscf_point = nscf_point.split('Writing')[0]
            qe_nscf.append(gen_lst(nscf_point.replace('\n', ''), ' ', float))
        if task == 'nscf':
            N_k, N_bands = np.array(qe_nscf).shape
            qe_nscf = np.reshape(qe_nscf[:N_k], (int(np.sqrt(N_k)), int(np.sqrt(N_k)), N_bands))
        else:
            qe_nscf = np.transpose(qe_nscf)
        return {'qe_'+task: qe_nscf}

    @change_directory('data_directory')
    def plot_bands(self, name='QE bands', bs=None, **kwargs):
        if bs == None:
            bs = self.DFT.gen_bs()
            kwargs['reset_color'] = True
        bands = np.load(self.prefix+'/qe_bands.npy')
        bs.add_bands(bands, name=name, **kwargs)
        return bs