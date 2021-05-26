from .decorators import change_directory, save
import numpy as np
import re

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_nscf(self, task='nscf'):
    nscf_file = read_file(self.prefix+'.'+task+'.out')
    nscf_data = nscf_file.split('End of band structure calculation')[1]
    nscf_data = gen_lst(nscf_data, 'k =')
    qe_nscf = []
    for i, line in enumerate(nscf_data):
        nscf_point = line.split('bands (ev):\n\n')[1].split('occupation')[0]
        nscf_point = nscf_point.split('Writing')[0]
        qe_nscf.append(gen_lst(nscf_point.replace('\n', ''), ' ', float))
    if task == 'nscf':
        N_k, N_bands = np.array(qe_nscf).shape
        qe_nscf = np.reshape(qe_nscf[:N_k], (int(np.sqrt(N_k)), int(np.sqrt(N_k)), N_bands))
    return {'qe_'+task: qe_nscf}

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_w90(self):
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
    ticks = [sum(self.spacing[:i]) for i in range(len(self.spacing)+1)]
    for i in range(len(ticks)-1):
        Δx = self.k_spacing[ticks[i+1]] - self.k_spacing[ticks[i]]
        Δn = n_ticks[i+1] - n_ticks[i]
        x_arr += [n*Δx/Δn+self.k_spacing[ticks[i]] for n in range(Δn)]
    x_arr += [self.k_spacing[-1]]
    hopping_file = read_file(self.prefix+'_hr.dat')
    hop_list = [gen_lst(lst, '   ', float) for lst in gen_lst(hopping_file, '\n')[8:-1]]

    return dict(w90_bands=w90_bands, w90_band_ticks=x_arr, w90_hopping=hop_list)

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_orbitals(self, orbital_dir=''):
    iso_list = []
    orbital_data = {}
    if orbital_dir != '':
        orbital_dir += '/'
    for i in range(1, self.w90_dic['num_wann']+1):
        f = read_file(self.prefix+'_'+str(i).zfill(5)+'.xsf')
        data = f.split('BEGIN_DATAGRID_3D_UNKNOWN\n')[1]
        data = data.split('\n')
        orbital_data[orbital_dir+'w90_N_grid'] = np.array([gen_lst(data[0], ' ', int)], dtype=int)
        orbital_data[orbital_dir+'wan_origins'] = gen_lst(data[1], ' ', float)
        orbital_data[orbital_dir+'w90_vec_span'] = [gen_lst(data[j], ' ', float) for j in range(2, 5)]
        iso_data = np.array([gen_lst(row, ' ', float) for row in data[5:-3]]).flatten()
        iso_data = iso_data.reshape(*np.flip(orbital_data[orbital_dir+'w90_N_grid']))
        iso_list.append(np.swapaxes(iso_data, 0, 2))
    orbital_data[orbital_dir+'w90_orbitals'] = np.array(iso_list)
    return orbital_data

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_relax(self):
    f = read_file(self.prefix+'relax.out')
    position_data = f.split('ATOMIC_POSITIONS (alat)\n')[1:]
    position_data = [item.split('\nEnd final')[0] for item in position_data]
    data_dic = {}
    for i, part in position_data:
        for row in [gen_lst(row, ' ') for row in part.split('\n')]:
            dic_key = 'iter_' + str(i) + '_' + row[0]
            if dic_key in data_dic:
                data_dic[dic_key].append(row[1:])
            else:
                data_dic[dic_key] = [row[1:]]
    return data_dic

@change_directory('data_directory')
@save
@change_directory('work_directory')
def extract_projwfc(self, proj_dir=''):
    f = read_file(proj_dir+self.prefix+'.projwfc.out')
    projwfc_data = f.split(' k =   ')
    states = gen_lst(projwfc_data[0], '\n     state #', separate_state, True)
    occupations = []
    bands = []
    k_points = []
    for k_point in projwfc_data[1:]:
        band_k = []
        band_data = k_point.split('\n    |psi|^2')[:-1]
        occupation = np.zeros((len(band_data), len(states)))
        for i, band in enumerate(band_data):
            ε, ψ = band.split(' eV ==== \n     psi = ')
            band_k.append(scrub_str(ε.split(') = ')[1]))
            ψ = gen_lst(ψ, '+', lambda x: scrub_str(x, '*'))
            for φ in ψ:
                occupation[i, int(φ[1])-1] = φ[0]
        occupations.append(occupation)
        bands.append(band_k)
        k_points.append(k_point)
    return {
        'projwfc/k_points': k_points, 
        'projwfc/bands': bands, 
        'projwfc/states': states, 
        'projwfc/occupations': occupations
    }

def scrub_str(string, char=None):
    if char == None:
        return float(re.sub("[^0-9.-]", "", string))
    else:
        return [float(re.sub("[^0-9.-]", "", x)) for x in string.split(char)]

def gen_lst(lst, str, func=lambda x: x, ignore_first=False):
    new_lst = []
    for i, item in enumerate(lst.split(str)):
        if (not empty(item)) and (i!=0 or not ignore_first):
            new_lst.append(func(item))
    return new_lst

def read_file(file_name):
    f = open(file_name, "r")
    file_content = f.read()
    f.close()
    return file_content

def separate_state(state):
    atom, q_num = state.split(', wfc')
    q_num = q_num.replace('= ', '=')
    atom = int(atom.split('atom')[1].split('(')[0])
    q_num = q_num.split('(')[1].split(')')[0]
    q_num = [scrub_str(q) for q in q_num.split(' ')]
    return [atom] + q_num

def empty(item):
    return item==[] or re.sub(r'[ \n]', '', str(item))==''