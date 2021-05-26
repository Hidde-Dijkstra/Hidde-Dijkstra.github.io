from .decorators import check

@check('qe')
def import_lattice(self, lattice):
    self.lattice = lattice
    self.qe_dic['&SYSTEM']['ibrav'] = 4
    self.qe_dic['&SYSTEM']['a'] = lattice.a
    self.qe_dic['&SYSTEM']['c'] = self.Δz
    self.qe_dic['&SYSTEM']['nat'] = len(lattice.unit_cell)
    self.qe_dic['&SYSTEM']['ntyp'] = len(lattice.atom_types)
    self.atoms = {}
    for atom in lattice.unit_cell:
        if atom.name not in self.atoms:
            self.atoms[atom.name] = atom.__dict__
            self.atoms[atom.name]['loc'] = [atom.position*lattice.a]
        else:
            self.atoms[atom.name]['loc'].append(atom.position*lattice.a)