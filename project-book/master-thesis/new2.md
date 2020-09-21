---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.9'
    jupytext_version: 1.5.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Monolayer WSe$_2$

+++

We base our tight binding model of monolayer WSe$_2$ (mWSe$_2$) on the three band model (TBM) by Liu et al. {cite}`three_band`. In the TBM of mWSe2 the hopping is modeled using only the tungsten sites, forming a triangular lattice in the $xy$ plane with real space unit vectors:
```{math}
\vec a_1=\begin{pmatrix}1\\ 0\end{pmatrix}\text{ and } \vec a_2=\frac12\begin{pmatrix}1\\ \sqrt{3}\end{pmatrix},    
```
and reciprocal lattice vectors:
```{math}
\vec b_1=2\pi\begin{pmatrix}1\\ -1/\sqrt{3}\end{pmatrix}\text{ and } \vec b_2=2\pi\begin{pmatrix}0\\ 2/\sqrt{3}\end{pmatrix},    
```
where length is in units of the lattice spacing $a$ between two tungsten atoms.

We consider one atom as the origin and label every other atom using the integers $i$ and $j$ which connect this atom to their origin with the vector $i
\vec a_1 + j\vec a_2$:

+++

```{glue:figure} fig:lattice
```

The symmetric point group of this triangular lattice is $D_3$ with elements $\{E, C_3, C^2_3,\sigma_v, \sigma_v', \sigma_v''\}$. The rotations are in the $xy$ plane, while the reflections are with respect to the $yz$, $xz$ and $xy$ planes. $C_6$ is not a symmetry element of this lattice since this operation exchanges triangles containing Se$_2$ and empty triangles.

We restrict the orbitals of tungsten to the $d$ orbitals $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$, which dictate the low energy $K$ points. The other orbitals are either at a different energy scale or do not mix with these orbitals due to the symmetries of the lattice. Since $C_6$ is not a symmetry element, $C_2=C_6\prod C_3$ is not either.

+++

```{glue:figure} fig:d_z2-d_xy
```

+++

```{admonition} GGA coupling constants in eV
Coupling constants for WSe$_2$, subscripts 1, 2 and 3 refer to the $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$ orbitals respectively. These constants are for NN hopping left to right with respect to the $x$-axis. $\varepsilon_i$ are the on-site energies where $\varepsilon_1=\varepsilon_2$ due to symmetry. $\lambda_\text{SOC}$ is the spin orbit coupling. These are the GGA constants from {cite}`three_band`.
| $t_1$| $t_2$| $t_3$| $t_{12}$| $t_{13}$| $t_{23}$| $\varepsilon_1$| $\varepsilon_3$| $\lambda_{SOC}$|
|:----:|:----:|:----:|:-------:|:-------:|:-------:|:--------------:|:--------------:|:--------------:|
|{glue:}`var:t_1`|{glue:}`var:t_2`|{glue:}`var:t_3`|{glue:}`var:t_12`|{glue:}`var:t_13`|{glue:}`var:t_23`|{glue:}`var:ε_1`|{glue:}`var:ε_3`|{glue:}`var:λ_SOC`|
where all energies are measured in eV. 
```
and collect them in the matrix $\mathcal E$:
```{math}
    \mathcal E = \begin{pmatrix}
        t_1&t_{12}&t_{13}\\
        -t_{12}&t_2&t_{23}\\
        t_{13}&-t_{23}&t_3
    \end{pmatrix}.
```

+++

```{bibliography} references.bib
```

```{code-cell} ipython3
import new1
```

```{code-cell} ipython3

```
