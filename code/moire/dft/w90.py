from .decorators import change_directory, save
from . import io_processing as io
from moire.bandstructure import BandStructure
import numpy as np

class W90:

    def __init__(self, DFT, w90_dic, spinors=False, N_super=3):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix       
        self.a = DFT.lattice.a
        self.a_vec = [self.a*DFT.a_1, self.a*DFT.a_2, DFT.Δz*DFT.a_3]
        self.Δa = DFT.middle * sum(self.a_vec)/2
        self.w90_dic = dict(guiding_centres=True, write_hr=True, wannier_plot=False,
                            bands_plot=True, spinors=spinors)
        self.w90_dic.update(w90_dic)
        if 'nbnd' in DFT.QE.qe_dic['&SYSTEM']:
            self.w90_dic['num_bands'] = DFT.QE.qe_dic['&SYSTEM']['nbnd']
        if self.w90_dic['wannier_plot']:
            self.w90_dic['wannier_plot_supercell'] = f'{N_super}, {N_super}, 1'
        self.orbital_centers = []
        self.a_vec = [self.a*DFT.a_1, self.a*DFT.a_2, DFT.Δz*DFT.a_3]

    def set_window(self, win_min, win_max, num_iter, froz_min=None, froz_max=None):
        self.w90_dic['dis_win_min'] = win_min
        self.w90_dic['dis_win_max'] = win_max
        self.w90_dic['dis_num_iter'] = num_iter
        if froz_max != None:
            self.w90_dic['dis_froz_max'] = froz_max
        if froz_min != None:
            self.w90_dic['dis_froz_min'] = froz_min

    @change_directory('work_directory')
    def write(self, pw2wan=True):
        if pw2wan:
            self.DFT.QE.write_pw2wan(self.w90_dic['wannier_plot'])
        file = [f'{setting} = {stringify(value)}'
                for setting, value in self.w90_dic.items()]
        file.append('\nbegin unit_cell_cart')
        file.extend([ f'  {io.join_grid_point(a_i)}' for a_i in self.a_vec])
        file.append('end unit_cell_cart\n\nBegin projections')
        file.extend([f"  {name}:  {'; '.join(atom.projections)}"
                     for name, atom in self.DFT.lattice.atom_types.items()
                     if atom.projections != []])
        file.append('End projections\n\nBegin atoms_cart\nang')
        file.extend([(f'  {atom.name.ljust(2)}'
                      f'  {io.join_grid_point(atom.position*self.a+self.Δa)}')
                     for atom in self.DFT.lattice])
        file.append('End atoms_cart\n\nBegin Kpoint_Path')
        sym_points = self.DFT.path.sym_points
        file.extend([(f'  {k1.tag} {io.join_grid_point(k1.vector)}'
                      f'  {k2.tag} {io.join_grid_point(k2.vector)}')
                     for k1, k2 in zip(sym_points, sym_points[1:])])
        file.append((f'End Kpoint_Path\n\nmp_grid = '
                     f'{self.DFT.N_grid}, {self.DFT.N_grid}, 1'))
        file.append((f'\nBegin kpoints\n{io.join_grid(self.DFT.k_grid, weight=False)}'
                     f'\nEnd kpoints\n'))
        io.write_file(f'{self.prefix}.win', '\n'.join(file))

    @change_directory('work_directory')
    def run(self, pp=True, n_cores=4):
        if pp:
            io.run(f'wannier90.x -pp {self.prefix}.win')
            self.DFT.QE.run('pw2wan', 'pw2wannier90.x', n_cores)
        io.run(f'wannier90.x {self.prefix}.win')

    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self):
        w90_data = io.delim(io.read_file(f'{self.prefix}.wout'),
                            'Final State',
                            'Sum of centres and spreads')
        band_file = io.read_file(f'{self.prefix}_band.dat')
        band_data = [[io.gen_lst(elem, ' ', float) 
                     for elem in io.gen_lst(lst, '\n')]
                    for lst in io.gen_lst(band_file, '\n  \n')]
        ticks_file = io.read_file(f'{self.prefix}_band.gnu')
        return dict(
            wan_centers = [io.gen_lst(io.delim(center, '(', ')'), ',', float) 
                           for center in io.gen_lst(w90_data, '\n')],
            wan_spreads = [float(center.split(' ')[-1])
                           for center in io.gen_lst(w90_data, '\n')],
            n_ticks = [np.where(np.array(band_data)[0, :, 0] >= k-10E-6)[0][0] 
                       for k in io.scrub_str(ticks_file.split('\n')[6], ',')],
            w90_bands = np.array(band_data)[:, :, 1]
        )
    
    @change_directory('data_directory')
    def plot_bands(self, name='W90 bands', bs=None, **kwargs):
        path = self.DFT.path.copy(N=np.load(f'{self.prefix}/n_ticks.npy'))
        if bs == None:
            bs = BandStructure(kwargs.get('cut'), path)
            kwargs['reset_color'] = True
        bands = np.load(f'{self.prefix}/w90_bands.npy')
        bs.add_bands(bands, name=name, path=path, **kwargs)
        return bs
   
def stringify(value):
    if type(value) == list:
        return ', '.join([str(item) for item in value])
    else:
        return str(value).lower()