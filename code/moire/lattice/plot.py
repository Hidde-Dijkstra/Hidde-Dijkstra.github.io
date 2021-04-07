from collections import defaultdict
import numpy as np
import plotly.graph_objects as go

class _DrawLattice:


    def __init__(self, lattice, W, H=None, D=None):
        if H == None:
            H = W
        if D == None:
            D = W
        self.W = W
        self.H = H
        self.D = D
        self.z0 = np.average([atom.position[2] for atom in lattice.unit_cell])
        self.unit_cell = lattice.unit_cell
        self.a_1 = lattice.a_1
        self.a_2 = lattice.a_2
        self.a = lattice.a

    def plot(self, plot_3d=False, atom_keys=[], bonds=False, show_z=False, **kwargs):
        general_dic = dict(
            mode='markers',
            showlegend=False,
            line_width=1
        )
        if 'all' in kwargs:
            for key, value in kwargs['all'].items():
                general_dic[key] = value
        fig = go.Figure()
        fig.update_xaxes(range=[-self.W/2, self.W/2])
        fig.update_yaxes(scaleanchor='x', scaleratio=1, range=[-1*self.H/2, self.H/2])
        if plot_3d:
            z_range = dict(range=[-1*self.D/2, self.D/2])
            fig.update_layout(scene=dict(zaxis=z_range))
        scatter_dic, position_dic = self._gen_plot_dics(general_dic, atom_keys)
        for atom_key, positions in position_dic.items():
            x, y, z = np.transpose(positions)
            scatter_dic[atom_key]['x'] = x
            scatter_dic[atom_key]['y'] = y
            if plot_3d:
                scatter_dic[atom_key]['z'] = z - self.z0
                fig.add_trace(go.Scatter3d(**scatter_dic[atom_key]))
            else:
                if show_z:
                    scatter_dic['marker']['color'] = z
                fig.add_trace(go.Scatter(**scatter_dic[atom_key]))
        return fig

    def _gen_plot_dics(self, general_dic, atom_keys):
        position_dic = defaultdict(list)
        scatter_dic = {}
        for atom in self.unit_cell:
            if atom_keys == [] or atom.key in atom_keys:
                if atom.key+atom.name not in scatter_dic:
                    scatter_dic[atom.key+atom.name] = {}
                    for key, value in general_dic.items():
                        scatter_dic[atom.key+atom.name][key] = value
                    scatter_dic[atom.key+atom.name]['marker'] = dict(
                        color=atom.color, 
                        size=atom.draw_radius,
                        line_width=general_dic['line_width']
                    )
                    scatter_dic[atom.key+atom.name]['name'] = atom.key
                for position in self._get_positions(atom):
                    position_dic[atom.key+atom.name].append(position)
        return scatter_dic, position_dic

    def _get_positions(self, atom):
        positions = []
        Δn = int(max(self.W, self.H)/np.linalg.norm(self.a_1)) + 2
        n_range = range(-Δn, Δn+1)
        for vec in [i*self.a_1+j*self.a_2 for i in n_range for j in n_range]:
            if self._in_grid(vec, atom):
                positions.append(atom.position+vec)
        return positions

    def _in_grid(self, vec, atom):
        origin = np.abs((atom.position+vec))
        return np.all((origin)[:2] < [self.W/2, self.H/2])  