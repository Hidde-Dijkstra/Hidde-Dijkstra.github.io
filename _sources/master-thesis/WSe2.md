---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.6.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Monolayer WSe$_2$


We base our tight binding model of monolayer WSe$_2$ (mWSe$_2$) on the three band model by Liu et al. {cite}`three_band`. In this model of mWSe2 we consider only the tungsten sites, which form a triangular lattice in the $xy$ plane: 


```{glue:figure} fig:lattice
Triangular lattice of WSe$_2$, the blue atoms respresent the tungsten atoms while the red the (di)selenide.
```

with real space unit vectors:
```{math}
\vec a_1=\begin{pmatrix}1\\ 0\end{pmatrix}\text{ and } \vec a_2=\frac12\begin{pmatrix}1\\ \sqrt{3}\end{pmatrix},    
```
and reciprocal lattice vectors:
```{math}
\vec b_1=2\pi\begin{pmatrix}1\\ -1/\sqrt{3}\end{pmatrix}\text{ and } \vec b_2=2\pi\begin{pmatrix}0\\ 2/\sqrt{3}\end{pmatrix},    
```
where length is in units of the lattice spacing $a$ between two tungsten atoms.


## Symmetries

The symmetric point group of this triangular lattice is $D_3$ with generators $\{C_3, \sigma\}$. $C_3$ is a rotation of $2\pi/3$ in the $xy$ plane and $\sigma$ is a reflection with respect to the bisector of $\vec a_1$ and $\vec a_2$. $C_6$ is not a symmetry element of this lattice since this operation exchanges triangles containing Se$_2$ and empty triangles.

We restrict the orbitals of tungsten to the $d$ orbitals $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$, which dictate the low energy $K$ points. The other orbitals are either at a different energy scale or do not mix with these orbitals due to the symmetries of the lattice. The tight binding hopping matrix between these orbitals in the hopping direction $\vec a_1$ writes:

```{math}
    \mathcal H_{\vec{a}_1} = \begin{pmatrix}
        t_1&t_{12}&t_{13}\\
        -t_{12}&t_2&t_{23}\\
        t_{13}&-t_{23}&t_3
    \end{pmatrix},
```

in the basis $(d_{x^2-y^2}\, d_{xy}\,d_{z^2})$.


```{admonition} GGA coupling constants in eV
GGA coupling constants for WSe$_2$ in the hopping direction $\vec{a}_1$, subscripts 1, 2 and 3 refer to the $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$ orbitals respectively. $\varepsilon_i$ are the on-site energies where $\varepsilon_1=\varepsilon_2$ due to symmetry. $\lambda_\text{SOC}$ is the spin orbit coupling. {cite}`three_band`.
| $t_1$| $t_2$| $t_3$| $t_{12}$| $t_{13}$| $t_{23}$| $\varepsilon_1$| $\varepsilon_3$| $\lambda_{SOC}$|
|:----:|:----:|:----:|:-------:|:-------:|:-------:|:--------------:|:--------------:|:--------------:|
|{glue:}`var:t_1`|{glue:}`var:t_2`|{glue:}`var:t_3`|{glue:}`var:t_12`|{glue:}`var:t_13`|{glue:}`var:t_23`|{glue:}`var:ε_1`|{glue:}`var:ε_3`|{glue:}`var:λ_SOC`|
where all energies are measured in eV. 
```



Surprisingly, $\mathcal H_{\vec a_1}$ is not symmetric. We recognize the lack of $C_6$ symmetry as the source of this irregularity. This is more easily understood when we consider the interaction between orbitals $d_{xy}$ and $d_{z^2}$:

```{glue:figure} fig:d_z2-d_xy
Representation of atomic orbitals $d_{xy}$ and $d_{z^2}$, the colors represent the different signs of the lobes.
```
Mirroring the hopping direction is equivalent to applying a reflection with respect to the $y$-axis. The $d_{xy}$ orbital changes sign under this operation so $t_{32}=-t_{23}$. If this reflection was a symmetry element of the lattice we would require $t_{32}=t_{23}$, forcing $t_{23}=0$. We can use a similar argument to exclude the $d_{xz}$ and $d_{yz}$ orbitals from our model, since for a monolayer we do have a reflection symmetry with respect to the $xy$ plane.


## Rotating orbitals

Having figured out the hopping terms for $\pm\vec a_1$, we need also formulate the coupling between orbitals at an angle with respect to the principal axis:

```{glue:figure} fig:rotate_orbital
Two $d_{x^2-y^2}$ orbitals at angle $\theta$ with respect to the principal axis of the orbitals.
```


```{bibliography} references.bib
```
