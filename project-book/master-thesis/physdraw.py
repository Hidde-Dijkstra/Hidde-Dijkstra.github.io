import drawSvg as draw
import numpy as np

class Index:
    
    def __init__(self, i, j):
        self.i = i
        self.j = j
        
    def __mul__(self, other):
        return self.i * other[0] + self.j * other[1]
    
    def __rmul__(self, other):
        return other * self
        
        
class Lattice:
    
    
    def __init__(self, a_1, a_2, W, H, θ=0):
        self.a = [a_1, a_2]
        self.W = W
        self.H = H
        self.unit_cell = []
        self.grid = []
        self.θ = θ
        
    def add_atom(self, atom):        
        N_1, N_2 = [1+min(int(self.W/np.abs(a[0]+0.00001)), int(self.H/np.abs(a[1]+0.00001))) for a in [rot_mat(self.θ)@a for a in self.a]]
        self.unit_cell.append(atom)
        self.grid.append([Index(i, j) for i in range(-N_1, N_1+1) for j in range(-N_2, N_2+1) if self.in_lattice(i, j, atom)])
        
    def in_lattice(self, i, j, atom):
        origin = np.abs(rot_mat(self.θ) @ (position+Index(i, j)*self.a))
        return np.all(origin < [self.W/2, self.H/2])
    
    def draw_lattice(self):
        container = draw.Drawing(self.W, self.H, origin="center", displayInline=False)
        for i, atom in enumerate(self.unit_cell):
            for grid_point in self.grid[i]:
                component_list = atom.draw_bonds(grid_point*self.a, self.θ)
                [container.append(component) for component in component_list]
        for i, atom in enumerate(self.unit_cell):
            if atom.atom_radius == None:
                atom.atom_radius = min([np.linalg.norm(a) for a in self.a]) / np.sqrt(len(self.unit_cell)) / 5 
            for grid_point in self.grid[i]:
                component_list = atom.draw_atom(grid_point*self.a, self.θ)
                [container.append(component) for component in component_list]
        return container
                
    def NN(self, atom_1, atom_2, bond_list):
        for bond in bond_list:
            atom_1.bonds.append(atom_2.position+Index(*bond)*self.a)
            if atom_1 != atom_2:
                atom_2.bonds.append(-atom_1.bonds[-1])
    


class LatticeAtom:
    
    
    def __init__(self, position_in_unit_cell, name=None, atom_color='blue', atom_radius=None):       
        self.position = np.array(position_in_unit_cell)
        self.name = name
        self.atom_color = atom_color
        self.atom_radius = atom_radius
        self.bonds = []
    
    def draw_bonds(self, displacement, θ, **kwargs):
        bond_components = []
        origin = rot_mat(θ) @ (displacement + self.position)
        for bond in self.bonds:
            bond_components.append(draw.Line(*origin, *(origin+bond), stroke='white', stroke_width=0.01, **kwargs))
        return bond_components
        
    def draw_atom(self, displacement, θ, **kwargs):
        atom_components = []
        origin = rot_mat(θ) @ (displacement + self.position)
        gradient = draw.RadialGradient(*origin, self.atom_radius)
        gradient.addStop(0, 'white', 1)
        gradient.addStop(1, self.atom_color, 1)
        atom_components.append(draw.Circle(*origin, self.atom_radius, fill=gradient, **kwargs))
        if self.name != None:
            atom_components.append(draw.Text(self.name, self.atom_radius, *origin, text_anchor='middle', alignment_baseline="central"))
        return atom_components 
        
class Orbital:  
    
    def lobe(self, color, rotate=0, translate=(0, 0), stroke="black", **kwargs):
        gradient = draw.RadialGradient(0, 100, 50)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(173.21, color, 0.7)
        transform = "translate(" + " ".join([str(i) for i in translate]) + ")\nrotate(" + str(rotate) + " 0 0)"
        return draw.Path(d="M 0,0 C -173.21,-200 173.21,-200 0,0 z", stroke=stroke, fill=gradient, transform=transform, **kwargs)
    
    def circle(self, color, ellipse=False, rotate=0, translate=(0, 0), stroke="black", **kwargs):
        gradient = draw.RadialGradient(0, 0, 50)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(173.21, color, 0.7)
        transform = "rotate(" + str(rotate) + " 0 0)\ntranslate(" + " ".join([str(i) for i in translate]) + ")"
        if ellipse:
            clip = draw.ClipPath()
            clip.append(draw.Ellipse(0, 0, 50, 12.5, transform=transform))
            return draw.Ellipse(0, 0, 100, 25, stroke=stroke, fill=gradient, transform=transform, **kwargs) 
        else:
            return draw.Circle(0, 0, 50, stroke=stroke, fill=gradient, transform=transform, **kwargs)

class d_xy(Orbital):
    
    def __init__(self, container, translate=(0, 0), rotate=0):
        container.append(self.lobe("dodgerblue", rotate=85+rotate, translate=translate))
        container.append(self.lobe("red", rotate=95+rotate, translate=translate))
        container.append(self.lobe("red", rotate=275+rotate, translate=translate))
        container.append(self.lobe("dodgerblue", rotate=265+rotate, translate=translate))
        
class d_z2(Orbital):
    
    def __init__(self, container, translate=(0, 0), rotate=0):
        container.append(self.lobe("dodgerblue", rotate=180+rotate, translate=translate))
        container.append(self.circle("red", ellipse=True, rotate=rotate, translate=translate))
        container.append(self.lobe("dodgerblue", rotate=rotate, translate=translate))

def rot_mat(θ):
    return np.array([[np.cos(θ), -np.sin(θ)], [np.sin(θ), np.cos(θ)]])