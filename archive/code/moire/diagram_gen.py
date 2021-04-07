# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.8.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Diagram drawing

import numpy as np
import plotly.graph_objects as go
import drawSvg as draw


# ## drawSvg
#
# The diagrams we use are based around the python module drawSvg which allows construction of svg diagrams straight from python. The main container for an image is always generated in the chapter file which is then filled with svg element groups which can pass on any arguments to its child elements (such as circles, lines etc.).

# ## Utility functions
#
# ### Rotation matrix
#
# A function which returns a rotation matrix in the xy plane in two or three dimensions.
#
# :::{dropdown} `rot_mat(θ, dim=2)`
# * `θ` rotation angle.
# * `dim` dimension of desired matrix, either 2 or 3.
# :::

def rot_mat(θ, dim=2):
    if dim == 2:
        return np.array([[np.cos(θ), -np.sin(θ)], [np.sin(θ), np.cos(θ)]])
    else:
        return np.array([[np.cos(θ), -np.sin(θ), 0], [np.sin(θ), np.cos(θ), 0], [0, 0, 1]])


# ### Arrow
#
# Draws an arrow between two points.
#
# :::{dropdown} `arrow(start, end, stroke_width=0.1, stroke='black', **kwargs)`
# * `start` starting point of the arrow.
# * `end` end point of the arrow with arrowhead.
# * `stroke_width` width of the arrow line.
# * `stroke` color of the arrow.
# :::

def arrow(start, end, stroke_width=0.1, stroke='black', **kwargs):
    start, end = np.array(start), np.array(end)
    Δx = 3
    my_arrow = draw.Marker(-1+Δx/4, -0.5, Δx/4, 0.5, scale=4, orient='auto')
    my_arrow.append(draw.Lines(-1+Δx/4, -0.5, -1+Δx/4, 0.5, Δx/4, 0, close=True, fill=stroke))
    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none',
              marker_end=my_arrow, **kwargs)
    t = 1 - stroke_width*Δx/np.linalg.norm(end-start)
    return p.M(*start).L(*(t*(end-start)+start))


# ### Index
#
# Object which functions as the index of an atomic site. We overload the multiplication operator to allow for shorthand of $ia_1+ja_2$. 
#
# ::::{dropdown} `Index`
# * :::{dropdown} `__init__(i, j)`
#     * `i` first index.
#     * `j` second index.
#     :::
# * :::{dropdown} `__mul__(other)`
#     * `other` a pair of two objects which we want to add with some integer weight, usually [$a_1$, $a_2$].
#     :::
# ::::

class Index:
 

    def __init__(self, i, j):
        self.i = i
        self.j = j
        
    def __mul__(self, other):
        return self.i * other[0] + self.j * other[1]
    
    def __rmul__(self, other):
        return self * other


# ## LatticeAtom
#
# We introduce a class which keeps track of the properties of individual atoms in the unit cell with methods:
#
# ::::{dropdown} `LatticeAtom`
# * :::{dropdown} `__init__(position_in_unit_cell, name=None, atom_color='blue', atom_radius=None)`
#     * `position_in_unit_cell` the location of the atom in the unit cell.
#     * `name` letter e.g. C for carbon indicating the type of atom, if None no name included.
#     * `atom_color` color to draw atom in.
#     * `atom_radius` size of atom to plot, if None no bonds are plotted.
#     :::
# * :::{dropdown} `draw_bonds(displacement, θ, **kwargs)`
#     draws atomic bonds around the atom.
#     * `displacement` translation with which to draw the bonds, usually a multiple of the lattice vectors.
#     * `θ` angle with which the lattice is rotated.
#     :::
# * :::{dropdown} `draw_atom(displacement, θ, **kwargs)`
#     draws the atom.
#     * `displacement` translation with which to draw the atom, usually a multiple of the lattice vectors.
#     * `θ` angle with which the lattice is rotated.
#     :::
# ::::

class LatticeAtom:
    
    
    def __init__(self, position_in_unit_cell, name=None, atom_color='blue', atom_radius=None):       
        self.position = np.array(position_in_unit_cell)
        self.name = name
        self.atom_color = atom_color
        self.atom_radius = atom_radius
        self.bonds = []
        self.bond_style = [] #
    
    def draw_bonds(self, displacement, θ, **kwargs):
        group = draw.Group()
        origin = rot_mat(θ) @ (displacement + self.position)
        for bond in self.bonds:
            destination = rot_mat(θ) @ (displacement+bond)
            group.append(draw.Line(*origin, *destination, stroke='black', stroke_width=0.01, **kwargs))
        return group
    
    def draw_atom(self, displacement, θ, **kwargs):
        group = draw.Group()
        origin = rot_mat(θ) @ (displacement + self.position)
        gradient = draw.RadialGradient(*origin, self.atom_radius)
        gradient.addStop(0, 'white', 1)
        gradient.addStop(1, self.atom_color, 1)
        group.append(draw.Circle(*origin, self.atom_radius, stroke='black', stroke_width=0.01, fill=gradient, **kwargs))
        if self.name != None:
            group.append(draw.Text(self.name, self.atom_radius, *origin, text_anchor='middle', alignment_baseline="central"))
        return group


# ## Lattice
#
# We introduce a class which collects all LatticeAtom objects in the unit cell and expands the lattice to fit into a diagram of size $W\times H$. 
#
# ::::{dropdown} `Lattice`
# * :::{dropdown} `__init__(a_1, a_2, W, H, θ=0)`
#     * `a_1` first lattice vector of 2d lattice.
#     * `a_2` second lattice vector of 2d lattice.
#     * `W` desired width of diagram.
#     * `H` desired height of diagram.
#     * `θ` rotation angle of lattice.
#     :::
# * :::{dropdown} `add_atom(atom)`
#     * `atom` LatticeAtom object to be added to the unit cell.
#     :::
# * :::{dropdown} `in_lattice(i, j, atom)`
#     checks if the position of `atom` plus $ia_1+ja_2$ lies within desired frame of the diagram.
#     * `i` first index.
#     * `j` second index.
#     * `atom` `LatticeAtom` object to be tested.
#     :::
# * :::{dropdown} `NN(atom_1, atom_2, bond_list, **kwargs)`
#     designate two ```{python}LatticeAtom``` objects as nearest neighbors
#     * `atom_1` first `LatticeAtom` object.
#     * `atom_2` second `LatticeAtom` object.
#     :::
# * :::{dropdown} `draw_lattice()`
#     returns svg group with all atoms and bonds in diagram frame. 
#     :::
# * :::{dropdown} `draw_lattice_vectors(vec_symbols=['a₁', 'a₂'], origin=(0, 0), centralize=True, stroke_width=0.1, color='black', **kwargs)`
#     returns svg group containing the two lattice vectors.
#     * `vec_symbols` list of two labels for the lattice vectors.
#     * `origin` coordinate from where to draw the lattice vectors.
#     * `centralize` if `True` translate the origin by $(a_1+a_2)/2$.
#     * `vec_symbols` list of two labels for the lattice vectors.
#     * `vec_symbols` list of two labels for the lattice vectors.
#     :::
# * :::{dropdown} `draw_unit_cell(origin=(0, 0), **kwargs)`
#     returns svg group containing dashed lines showing the unit cells.
#     * `origin` coordinate from where to start the unit cell.
#     ::: 
# ::::

class Lattice:
    
    unit_cell_NN = [Index(i, j) for i in [-1, 0, 1] for j in [-1, 0, 1]]
    
    def __init__(self, a_1, a_2, W, H, θ=0):
        self.a = [np.array(a_1), np.array(a_2)]
        self.dim = len(a_1)
        self.W = W
        self.H = H
        self.unit_cell = []
        self.grid = []
        self.θ = -θ
        
    def add_atom(self, atom):        
        N_1, N_2 = [1+min(int(self.W/np.abs(a[0]+0.00001)), int(self.H/np.abs(a[1]+0.00001))) for a in [rot_mat(self.θ, self.dim)@a for a in self.a]]
        self.unit_cell.append(atom)
        if atom.atom_radius == None:
            atom.atom_radius = min([np.linalg.norm(a) for a in self.a]) / min(N_1, N_2) / 5 
        self.grid.append([Index(i, j) for i in range(-N_1, N_1+1) for j in range(-N_2, N_2+1) if self.in_lattice(i, j, atom)])
        
    def in_lattice(self, i, j, atom):
        origin = np.abs(rot_mat(self.θ, self.dim) @ (atom.position+Index(i, j)*self.a))
        return np.all((origin-atom.atom_radius)[:2] < [self.W/2, self.H/2])
        
    def NN(self, atom_1, atom_2, bond_list, **kwargs):
        #atom_1.bond_style.append(kwargs)
        for bond in bond_list:
            atom_1.bonds.append(atom_2.position+Index(*bond)*self.a)
            if atom_1 != atom_2:
                #atom_2.bond_style.append(kwargs)
                #atom_2.bonds.append(-atom_1.bonds[-1])
                pass
                
    def draw_lattice(self):
        group = draw.Group()
        for i, atom in enumerate(self.unit_cell):
            for grid_point in self.grid[i]:
                group.append(atom.draw_bonds(grid_point*self.a, self.θ), z=0)                
                group.append(atom.draw_atom(grid_point*self.a, self.θ), z=1)
        return group
               
    def draw_lattice_vectors(self, vec_symbols=['a₁', 'a₂'], origin=(0, 0), centralize=True, stroke_width=0.1, color='black',
                             **kwargs):
        rot = rot_mat(self.θ)
        group = draw.Group()
        if centralize:
            origin += sum(self.a) / 2
        group.append(arrow(rot@origin, rot@(origin+self.a[0]), stroke_width=stroke_width, stroke=color, **kwargs))
        group.append(arrow(rot@origin, rot@(origin+self.a[1]), stroke_width=stroke_width, stroke=color, **kwargs))
        group.append(draw.Text(vec_symbols[0], stroke_width*10, *(rot@(origin+self.a[0])), fill=color))
        group.append(draw.Text(vec_symbols[1], stroke_width*10, *(rot@(origin+self.a[1])), fill=color))
        return group
    
    def draw_unit_cell(self, origin=(0, 0), **kwargs):
        rot = rot_mat(self.θ)
        group = draw.Group()
        for i in range(2):
            N = int(np.sqrt(self.H**2+self.W**2)/np.linalg.norm(self.a[1-i])) + 1
            for j in range(-N, N+1):
                vector = np.sqrt(self.H**2+self.W**2) * self.a[i]/np.linalg.norm(self.a[i])
                group.append(draw.Line(*(rot@(origin-vector+j*self.a[1-i])), *(rot@(origin+vector+j*self.a[1-i])), **kwargs))
        return group


# ## Lattice3d

# A child class of `Lattice` which supports 3d plotting to generate a plotly lattice.
# ::::{dropdown} `Lattice3d`
# * :::{dropdown} `gen_NN(R)`
#     generates a list of nearest neighbors for each atom in the unit cell within some radius $R$.
#     * `R` radius within which an atom is considered a nearest neighbor.
#     :::
# * :::{dropdown} `draw_3d(fig, origin=(0, 0, 0))`
#     adds a lattice to a plotly figure.
#     * `fig` 3d plotly figure to which we wish to add a lattice.
#     * `origin` origin of the lattice.
#     :::
# ::::

class Lattice3d(Lattice):
    
    def gen_NN(self, R):
        try:
            len(R)
        except:
            R = np.ones(len(self.unit_cell))*R
        for i in range(len(self.unit_cell)):
            for j in range(i, len(self.unit_cell)):
                for index in self.unit_cell_NN:
                    ΔR = self.unit_cell[j].position+self.a*index-self.unit_cell[i].position
                    if 10**-6 < np.linalg.norm(ΔR) < R[i]:
                        self.unit_cell[i].bonds.append(self.unit_cell[j].position+self.a*index)
                        if i != j:
                            self.unit_cell[j].bonds.append(self.unit_cell[i].position-self.a*index)
            
    def draw_3d(self, fig, origin=(0, 0, 0)):
        for i, atom in enumerate(self.unit_cell):
            marker=dict(color=atom.atom_color)
            x, y, z = np.array([self.a*index+atom.position+origin for index in self.grid[i]]).T
            for bond in atom.bonds:
                for j in range(len(x)):
                    x_2, y_2, z_2 = self.a*self.grid[i][j]+bond+origin
                    fig.add_trace(go.Scatter3d(x=[x[j], x_2], y=[y[j], y_2], z=[z[j], z_2], 
                                           mode='lines', marker=dict(color='black'), showlegend=False,
                                              name=atom.name))
            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=marker, showlegend=False, name=''))


# ## Orbital
#
# A class to draw orbital diagrams.
# ::::{dropdown} `Orbital`
# * :::{dropdown} `lobe(color, rotate=0, translate=(0, 0), stroke="black", **kwargs)`
#     a base lobe with which to construct more complex orbitals.
#     * `color`: color of the lobe.
#     * `rotate`: angle with which to rotate the lobe.
#     * `translate`: vector with which to translate the lobe.
#     * `stroke`: color of the border of lobe.
#     :::
# * :::{dropdown} `circle(color, rotate=0, translate=(0, 0), stroke="black", ellipse=False, **kwargs)`
#     a base circle with which to construct more complex orbitals.
#     * `color` color of the circle.
#     * `rotate` angle with which to rotate the circle.
#     * `translate` vector with which to translate the circle.
#     * `stroke` color of the border of circle.
#     * `ellipse` if * `True` the circle is deformed to an ellipse.
#     :::
# * :::{dropdown} `d_xy(translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red", **kwargs)`
#     * `translate` vector with which to translate the orbital.
#     * `rotate` angle with which to rotate the orbital.
#     * `neg_color` color of the negative parts of the orbital.    
#     * `pos_color` color of the positive parts of the orbital.
#     :::
# * :::{dropdown} `d_z2(translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red", **kwargs)`
#     * `translate` vector with which to translate the orbital.
#     * `rotate` angle with which to rotate the orbital.
#     * `neg_color` color of the negative parts of the orbital.    
#     * `pos_color` color of the positive parts of the orbital.
#     :::
# * :::{dropdown} `d_x2y2(translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red", **kwargs)`
#     * `translate` vector with which to translate the orbital.
#     * `rotate` angle with which to rotate the orbital.
#     * `neg_color` color of the negative parts of the orbital.    
#     * `pos_color` color of the positive parts of the orbital.
#     :::
# ::::

class Orbital:  
    
    
    def lobe(self, color, rotate=0, translate=(0, 0), stroke="black", **kwargs):
        gradient = draw.RadialGradient(0, 1, 0.5)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(np.sqrt(3), color, 0.7)
        transform = "translate(" + " ".join([str(i) for i in translate]) + ")\nrotate(" + str(rotate) + " 0 0)"
        my_path = "M 0,0 C " + str(-np.sqrt(3)) + ",-2 " + str(np.sqrt(3)) +",-2 0,0 z"
        return draw.Path(d=my_path, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs)
    
    def circle(self, color, rotate=0, translate=(0, 0), stroke="black", ellipse=False, **kwargs):
        gradient = draw.RadialGradient(0, 0, 0.5)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(np.sqrt(3), color, 0.7)
        transform = "rotate(" + str(rotate) + " 0 0)\ntranslate(" + " ".join([str(i) for i in translate]) + ")"
        if ellipse:
            clip = draw.ClipPath()
            clip.append(draw.Ellipse(0, 0, 0.5, 0.125, transform=transform))
            return draw.Ellipse(0, 0, 1, 0.25, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs) 
        else:
            return draw.Circle(0, 0, 0.5, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs)
    
    def d_xy(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
             **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=85+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=95+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=275+rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=265+rotate, translate=translate))
        return group
    
    def d_z2(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
             **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=180+rotate, translate=translate))
        group.append(self.circle(pos_color, ellipse=True, rotate=rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=rotate, translate=translate))
        return group
    
    def d_x2y2(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
               **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=180+rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=90+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=270+rotate, translate=translate))
        return group


