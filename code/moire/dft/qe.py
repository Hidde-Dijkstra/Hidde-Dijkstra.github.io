from dataclasses import dataclass
from .decorators import change_directory, save
from . import io_processing as io
from moire.bandstructure import BandStructure
import numpy as np

class QE:

    def __init__(self, DFT, qe_dic):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.a = DFT.lattice.a
        self.a_vec = [self.a*DFT.a_1, self.a*DFT.a_2, DFT.Δz*DFT.a_3]
        self.Δa = DFT.middle * sum(self.a_vec)/2
        self.qe_dic = {'&CONTROL': {'prefix': DFT.prefix,
                                    'outdir': './out',
                                    'verbosity': 'high'},
                       '&SYSTEM': {'assume_isolated': '2D',
                                   'nat': len(DFT.lattice.unit_cell),
                                   'ntyp': len(DFT.lattice.atom_types),
                                   'ibrav': 4,
                                   'a': DFT.lattice.a,
                                   'c': DFT.Δz},
                       '&ELECTRONS': {},
                       '&IONS': {}}
        for key, dic in qe_dic.items():
            self.qe_dic[key].update(dic)
            
            
    @change_directory('work_directory')
    def write(self, qe_type, k_auto=9, spinors=False):
        self.qe_dic['&CONTROL']['calculation'] = qe_type
        ecutrho, ecutwfc, atom_species = self.pseudo(spinors)
        self.qe_dic['&SYSTEM'].update(dict(ecutwfc=ecutwfc,
                                           ecutrho=ecutrho,
                                           lspinorb=spinors, 
                                           noncolin=spinors))
        file = []
        for section_name, section in self.qe_dic.items():
            file.append(section_name)
            file.extend([f'  {setting} = {stringify(value)}'
                         for setting, value in section.items()])
            file.append('/')
        file.append('ATOMIC_SPECIES')
        file.extend(atom_species)
        file.append('ATOMIC_POSITIONS angstrom')
        file.extend([(f'  {atom.name.ljust(2)}'
                      f'  {io.join_grid_point(atom.position*self.a+self.Δa)}')
                     for atom in self.DFT.lattice])
        if qe_type == 'scf' or qe_type == 'relax':
            file.append(f'K_POINTS automatic\n  {k_auto} {k_auto} 1 1 1 1')
        elif qe_type == 'nscf':
            file.append(f'K_POINTS crystal\n{len(self.DFT.k_grid)}')
            file.append(io.join_grid(self.DFT.k_grid))
        elif qe_type == 'bands':
            file.append(f'K_POINTS crystal\n{len(self.DFT.path)}')
            file.append(io.join_grid([k.vector for k in self.DFT.path]))
        io.write_file(f'{self.prefix}.{qe_type}.in', '\n'.join(file))
        
    def pseudo(self, spinors):
        ecut = ECUT().update(self.qe_dic['&SYSTEM'].get('ecutwfc', 0),
                             self.qe_dic['&SYSTEM'].get('ecutrho', 0))
        mode = 'relativistic' if spinors else 'standard'
        atom_species = []
        for name, atom in self.DFT.lattice.atom_types.items():
            pseudo = atom.pseudo.get(mode, atom.pseudo.get('standard'))
            ecut.update(pseudo.get('ecutwfc', 0), 
                        pseudo.get('ecutrho', 0))
            atom_species.append((f"  {name.ljust(2)}"
                                 f"  {str(atom.weight).ljust(8, '0')}"
                                 f"  '{pseudo.get('file', '')}'"))
        return ecut.rho, ecut.wfc, atom_species

    @change_directory('work_directory')
    def write_pw2wan(self, wannier_plot):
        control_dic = self.qe_dic['&CONTROL']
        file = [f"&inputpp",
                f"  outdir = '{control_dic['outdir']}'",
                f"  prefix = '{self.prefix}'",
                f"  seedname = '{self.prefix}'",
                f"  write_mmn = .true.",
                f"  write_amn = .true.",
                f"  write_unk = .{str(wannier_plot).lower()}.",
                f"/\n"]
        io.write_file(f'{self.prefix}.pw2wan.in', '\n'.join(file))

    @change_directory('work_directory')
    def run(self, task, executor='pw.x', n_cores=4):
        io.run((f'mpirun -n {n_cores} {executor}'
                f' < {self.prefix}.{task}.in'
                f' > {self.prefix}.{task}.out'))
    
    @change_directory('data_directory')
    @save
    @change_directory('work_directory')
    def extract(self, task='nscf'):
        nscf_file = io.read_file(f'{self.prefix}.{task}.out')
        nscf_data = nscf_file.split('End of band structure calculation')[1]
        qe_nscf = [io.gen_lst(io.delim(line, 'bands (ev):', 'occupation', 
                                       'Writing').replace('\n', ''), 
                              ' ', float)
                   for line in io.gen_lst(nscf_data, 'k =')]
        if task == 'nscf':
            N_k, N_bands = np.array(qe_nscf).shape
            qe_nscf = np.reshape(qe_nscf[:N_k], (int(np.sqrt(N_k)), int(np.sqrt(N_k)), N_bands))
        else:
            qe_nscf = np.transpose(qe_nscf)
        return {f'qe_{task}': qe_nscf}

    @change_directory('data_directory')
    def plot_bands(self, name='QE bands', bs=None, N=None, ΔN=None, **kwargs):
        path = self.DFT.path.copy(N=N, ΔN=ΔN)
        if bs == None:
            bs = BandStructure(cut=kwargs.get('cut'), path=path)
            kwargs['reset_color'] = True
        bands = np.load(f'{self.prefix}/qe_bands.npy')
        bs.add_bands(bands, name=name, path=path, **kwargs)
        return bs
    
def stringify(value):
    if type(value) == str:
        return f"'{value}'"
    elif type(value) == bool:
        return str(value).lower()
    else:
        return str(value)

@dataclass
class ECUT:
    wfc = 0
    rho = 0
    
    def update(self, wfc, rho):
        self.wfc = wfc if wfc > self.wfc else self.wfc
        self.rho = rho if rho > self.rho else self.rho
        return self