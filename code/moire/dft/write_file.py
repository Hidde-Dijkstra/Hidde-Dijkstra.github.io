from .decorators import change_directory, check

@check('qe')
@change_directory('work_directory')
def write_qe(self, qe_type, k_auto=9):
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
        for a_i in self.a + [(0, 0, self.Δz)]:
            qe_file += '  ' + join_grid_point(a_i) + '\n'
    qe_file += 'ATOMIC_SPECIES\n'
    for atom_name, atom in self.atoms.items():
        qe_file += '  ' + '   '.join([
            atom_name.ljust(2), 
            str(atom['weight']).ljust(6, '0'), 
            "'"+atom['pseudo_potential']+"'\n"
        ]) 
    qe_file += 'ATOMIC_POSITIONS angstrom\n'
    for atom_name, atom in self.atoms.items():
        for site in atom['loc']:
            qe_file += '  ' + atom_name.ljust(2) + '   ' + join_grid_point(site) + '\n'
    if qe_type == 'scf' or qe_type == 'relax':
        qe_file += 'K_POINTS automatic\n  '+str(k_auto)+' '+str(k_auto)+' 1 1 1 1'
    elif qe_type == 'nscf':
        qe_file +=  'K_POINTS crystal\n' + str(len(self.k_grid)) + '\n'
        qe_file += join_grid(self.k_grid)
    elif qe_type == 'bands':
        qe_file +=  'K_POINTS crystal\n' + str(len(self.band_plot.k_path)) + '\n'
        qe_file += join_grid(self.band_plot.k_path)
    write_file(self.prefix+'.'+qe_type+'.in', qe_file)

@check('w90')
@change_directory('work_directory')
def write_w90(self, pw2wan=True):
    if pw2wan:
        write_pw2wan(self.qe_dic['&CONTROL'], self.w90_dic['wannier_plot'])
    w90_file = ''
    for setting, value in self.w90_dic.items():
        if type(value) == list:
            w90_file += setting + ' = ' + ', '.join([str(item) for item in value]) + '\n'
        else:
            w90_file += setting + ' = ' + str(value).lower() + '\n'
    w90_file += '\nbegin unit_cell_cart\n'
    for a_i in self.a:
        w90_file += '  ' + join_grid_point(a_i) + '\n'
    w90_file += 'end unit_cell_cart\n\nBegin projections'
    for atom_name, atom in self.atoms.items():
        if 'projections' in atom:
            w90_file += '\n  ' + atom_name +':  ' + '; '.join(atom['projections'])
    w90_file += '\nEnd projections\n\nBegin atoms_cart\nang\n'
    for atom_name , atom in self.atoms.items():
        for site in atom['loc']:
            w90_file += '  ' + atom_name + '   ' + join_grid_point(site) + '\n'
    w90_file += 'End atoms_cart\n\nBegin Kpoint_Path\n'
    for i in range(len(self.k_tags)-1):
        w90_file += '  ' + self.k_tags[i] + ' ' + join_grid_point(self.k_list[i]) + '   '
        w90_file += self.k_tags[i+1] + ' ' + join_grid_point(self.k_list[i+1]) + '\n'
    w90_file += 'End Kpoint_Path\n\nmp_grid = ' + str(self.N_grid) +', ' + str(self.N_grid) +', 1'
    w90_file += '\n\nBegin kpoints\n' + join_grid(self.k_grid, weight=False) + 'End kpoints'
    write_file(self.prefix+'.win', w90_file)

def write_pw2wan(control_dic, wannier_plot):
    pw2wan_file = "&inputpp\n   outdir = '" + control_dic['outdir'] +"'\n   "
    pw2wan_file += "prefix = '" + control_dic['prefix'] +"'\n   "
    pw2wan_file += "seedname = '" + control_dic['prefix'] +"'\n   "
    if wannier_plot:
        pw2wan_file += 'write_unk = .true.\n'
    pw2wan_file  += 'write_mmn = .true.\n   write_amn = .true.\n/\n'
    write_file(control_dic['prefix']+'.pw2wan.in', pw2wan_file)
    
def write_file(file_name, file_content):
    f = open(file_name, "w")
    f.write(file_content)
    f.close()

def join_grid_point(grid_point):
    return '   '.join(['{:.9f}'.format(item)[0:9] for item in grid_point])

def join_grid(k_grid, weight=True):
    grid_str = ''
    for grid_point in k_grid:
        grid_str += '  ' + join_grid_point(grid_point) + ' 1'*weight+'\n'
    return grid_str