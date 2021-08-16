from .decorators import change_directory, save
from . import io_processing as io
import numpy as np
from moire.plot.heatmap import gen_heatmap
import warnings

class RELAX:
    
    def __init__(self, DFT, relax_dir=''):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.dir = relax_dir
        self.a = DFT.QE.qe_dic['&SYSTEM']['a']
        
    def write(self, input_dft='vdw-df-cx', k_auto=1):
        self.DFT.QE.qe_dic['&SYSTEM']['input_dft'] = input_dft
        self.DFT.QE.write('relax', k_auto=k_auto)
        del self.DFT.QE.qe_dic['&SYSTEM']['input_dft']
        
    def run(self, n_cores=4):
        self.DFT.QE.run('relax', n_cores=n_cores)
    
    @change_directory('work_directory') 
    def extract(self, save=True):
        relax_file = io.read_file(f'{self.prefix}.relax.out')
        positions = [[io.gen_lst(line, ' ', float, ignore_first=True)
                      for line in split_relax(iteration)]
                     for iteration in relax_file.split('ATOMIC_POSITIONS')[1:]]
        atoms = [line.split(' ')[0] for line 
                 in split_relax(relax_file.split('ATOMIC_POSITIONS')[1])]
        if save:
            self._save_extract(positions, atoms)
        else:
            return positions, atoms
    
    @change_directory('data_directory')
    @save  
    def _save_extract(self, positions, atoms):
        return {f'{self.dir}/positions': positions,
                f'{self.dir}/atoms': np.array(atoms, dtype=str)}
    
    @change_directory('data_directory') 
    def plot(self, layers=[dict()], tol=0.1, iteration=-1, **kwargs):
        position_list = np.load(f'{self.prefix}/{self.dir}/positions.npy')
        atoms = np.load(f'{self.prefix}/{self.dir}/atoms.npy')
        tags = [layer.get('tag', f'layer {i+1}') 
                for i, layer in enumerate(layers)]
        positions = [[position
                      for atom, position in zip(atoms, position_list[iteration]) 
                      if (np.abs(position[2]-layer.get('z', position[2])) < tol 
                          and atom == layer.get('atom', atom))]
                     for layer in layers]
        return gen_heatmap(slider=dict(args=dict(positions=positions), tags=tags), 
                           a=self.a, **kwargs)

    def update_lattice(self, extract=False, tol=0.5):
        if extract:
            positions, atoms = self.extract(save=False)
        else:
            positions = np.load(f'{self.dir}/positions.npy')
            atoms = np.load(f'{self.dir}/atoms.npy')
        for atom, position, name in zip(self.DFT.lattice, 
                                        np.array(positions[-1])/self.a, 
                                        atoms):
            if (name == atom.name 
                and np.linalg.norm(atom.position-position) < tol):
                atom.position = position
            else:
                warnings.warn((f'Atom {atom.name} does not match {name}'
                               if atom.name != name else
                               f'Distance {atom.position} and {position} '
                               f'exceeds tolerance for relaxation.'))            
        
def split_relax(block):
    block = block.split('Writing output data file')[0]
    return io.gen_lst(block.split('End final coordinates')[0], '\n', 
                      ignore_first=True)