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


## Lattice
We base our tight binding model of monolayer WSe$_2$ (mWSe$_2$) on the three band model by Liu et al. {cite}`three_band`. In this model of mWSe2 we consider only the tungsten sites, which form a triangular lattice in the $xy$ plane: 


```{glue:figure} fig:lattice
Triangular lattice of WSe$_2$, the blue atoms respresent the tungsten atoms while the red the (di)selenide.
```

with real space unit vectors:

$$
\vec a_1=\begin{pmatrix}1\\ 0\end{pmatrix}\text{ and } \vec a_2=\frac12\begin{pmatrix}1\\ \sqrt{3}\end{pmatrix},    
$$

and reciprocal lattice vectors:

$$
\vec b_1=2\pi\begin{pmatrix}1\\ -1/\sqrt{3}\end{pmatrix}\text{ and } \vec b_2=2\pi\begin{pmatrix}0\\ 2/\sqrt{3}\end{pmatrix},    
$$

where length is in units of the lattice spacing $a$ between two tungsten atoms.


## Symmetries

The symmetric point group of this triangular lattice is $D_3$ with generators $\{C_3, \sigma\}$. $C_3$ is a rotation of $2\pi/3$ in the $xy$ plane and $\sigma$ is a reflection with respect to the bisector of $\vec a_1$ and $\vec a_2$. $C_6$ is not a symmetry element of this lattice since this operation exchanges triangles containing Se$_2$ and empty triangles.

We restrict the orbitals of tungsten to the $d$ orbitals $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$, which dictate the low energy $K$ points. The other orbitals are either at a different energy scale or do not mix with these orbitals due to the symmetries of the lattice. The tight binding hopping matrix between these orbitals in the hopping direction $\vec a_1$ writes:

$$
    \mathcal H_0 = \begin{pmatrix}
        t_1&t_{12}&t_{13}\\
        -t_{12}&t_2&t_{23}\\
        t_{13}&-t_{23}&t_3
    \end{pmatrix},
$$

in the basis $(d_{x^2-y^2}\, d_{xy}\,d_{z^2})$.


````{admonition} GGA coupling constants in eV
GGA coupling constants for WSe$_2$ in the hopping direction $\vec{a}_1$, subscripts 1, 2 and 3 refer to the $d_{x^2-y^2}$, $d_{xy}$ and $d_{z^2}$ orbitals respectively. $\varepsilon_i$ are the on-site energies where $\varepsilon_1=\varepsilon_2$ due to symmetry. $\lambda_\text{SOC}$ is the spin orbit coupling. {cite}`three_band`.
```{list-table}
:header-rows: 1
* - $t_1$
  - $t_2$
  - $t_3$
  - $t_{12}$
  - $t_{13}$
  - $t_{23}$
  - $\varepsilon_1$
  - $\varepsilon_3$
  - $\lambda_{SOC}$
* - {glue:}`var:t_1`
  - {glue:}`var:t_2`
  - {glue:}`var:t_3`
  - {glue:}`var:t_12`
  - {glue:}`var:t_13`
  - {glue:}`var:t_23`
  - {glue:}`var:ε_1`
  - {glue:}`var:ε_3`
  - {glue:}`var:λ_SOC`
```
where all energies are measured in eV. 
````


Surprisingly, $\mathcal H_0$ is not symmetric. We recognize the lack of $C_6$ symmetry as the source of this irregularity. This is more easily understood when we consider the interaction between orbitals $d_{xy}$ and $d_{z^2}$:

```{glue:figure} fig:d_z2-d_xy
Representation of atomic orbitals $d_{xy}$ and $d_{z^2}$, the colors represent the different signs of the lobes. We see that the $d_{z^2}$ orbital is symmetric under reflection while the $d_{xy}$ are antisymmetric.
```
Mirroring the hopping direction is equivalent to applying a reflection with respect to the $y$-axis. The $d_{xy}$ orbital changes sign under this operation so $t_{32}=-t_{23}$. If this reflection was a symmetry element of the lattice we would require $t_{32}=t_{23}$, forcing $t_{23}=0$. We can use a similar argument to exclude the $d_{xz}$ and $d_{yz}$ orbitals from our model, since for a monolayer we do have a reflection symmetry with respect to the $xy$ plane.


## Rotating orbitals

Having figured out the hopping terms for $\pm\vec a_1$, we need also formulate the coupling between orbitals at an angle with respect to the principal axis:

```{glue:figure} fig:rotate_orbital
Two $d_{x^2-y^2}$ orbitals at angle $\theta$ with respect to the principal axis of the orbitals.
```

We can however decompose any misaligned $d$ orbital as a linear sum of all $d$ orbitals aligned with the hopping axis. Since $d_{z^2}$ is symmetric under rotation in the $xy$ plane we can exclude this orbital from our analysis. We assume the existence of some matrix $R(\theta) = R'(\theta)\oplus 1$ which governs rotations in our basis and apply it to the $d_{x^2-y^2}$ orbital:

$$
\begin{align*}
R'(\theta)|{x^2-y^2}\rangle &= |{(x\cos\theta +y\sin\theta )^2-(y\cos\theta -x\sin\theta )^2}\rangle\\
&=|{x^2\cos^2\theta +y^2\sin^2\theta +xy\sin2\theta-y^2\cos^2\theta -x^2\sin^2\theta +xy\sin2\theta}\rangle\\
&=\cos2\theta|x^2-y^2\rangle + \sin2\theta|2xy\rangle,
\end{align*}
$$

which leads us to express $R'(\theta)$ in matrix form as:

$$
R'(\theta)=\begin{pmatrix}\cos2\theta&-\sin2\theta\\ \sin2\theta&\cos2\theta\end{pmatrix}\rightarrow \mathcal R=R'(\pi/3)\oplus 1=\begin{pmatrix}-1/2&-\sqrt{3}/2\\ \sqrt{3}/2&-1/2\end{pmatrix} \oplus 1.
$$

Powers of $\mathcal R$ now rotate the orbitals towards the respective hopping axes. Next we need also rotate $\mathcal H_0$ to adhere to the lattice symmetries. We use powers of $C_6$ which we express as $\text{Diag}([1, -1, 1])$ in accordance with the arguments of the previous part on lattice symmetry. 

We now give the hopping matrix for each hopping vector in: 

$$A = \{ \vec a_1, \vec a_2, -\vec a_1+\vec a_2, -\vec a_1, -\vec a_2, \vec a_1-\vec a_2 \} ,$$ 

where enumeration starts at zero:

$$
\mathcal H_i = \mathcal R^i C_6^i \mathcal H_0C_6^i \mathcal R^{-i}.
$$


## Hamiltonian

The total hopping Hamiltonian is now the sum of all $\mathcal H_i$ multiplied by the appropriate complex phases:

$$
\mathcal H_{hop}(\mathbf k) = \sum_{i=0}^5\mathcal H_i \text e^{i \mathbf k\cdot A_i} =    \begin{pmatrix}
        h_1&h_{12}&h_{13}\\ h_{12}^\dagger&h_2&h_{23}\\
        h_{13}^\dagger&h_{23}^\dagger&h_3
    \end{pmatrix} ,
$$

````{dropdown} with: (click dropdown for justification)

We start by noticing that all rotation matrices share a quite similar form:

```{list-table} The matrix expressions for $\mathcal R^iC_6^i$ excluding the $d_{z^2}$ component which is block diagonal.
:header-rows: 1
:name: rotation-matrices
* - $I_3$
  - $\mathcal RC_6$
  - $\mathcal R^2C_6^2=\mathcal R^2$
  - $\mathcal R^3C_6^3=C_6$
  - $\mathcal R^4C_6^4=\mathcal R$
  - $\mathcal R^5C_6^5=\mathcal R^2C_6$
* - $I_2$
  - $\begin{pmatrix}-1/2&\sqrt{3}/2\\\sqrt{3}/2&1/2 \end{pmatrix}$
  - $\begin{pmatrix}-1/2&\sqrt{3}/2\\-\sqrt{3}/2&-1/2 \end{pmatrix}$
  - $\begin{pmatrix}1&0\\0&-1 \end{pmatrix}$
  - $\begin{pmatrix}-1/2&-\sqrt{3}/2\\\sqrt{3}/2&-1/2 \end{pmatrix}$
  - $\begin{pmatrix}-1/2&-\sqrt{3}/2\\-\sqrt{3}/2&1/2 \end{pmatrix}$ 
```
As such we can express all of them as some variant of:

$$
\mathcal{\tilde R}(\alpha, \alpha', \beta, \beta')=\begin{pmatrix}\alpha&\beta&0\\\beta'&\alpha'&0\\0&0&1\end{pmatrix}
$$

where $\alpha^2=\alpha'^2$, $\beta^2=\beta'^2$ and $\alpha\beta'=-\alpha'\beta$. This rotation allows a general expression for the matrix components of $\mathcal H_i$ where the constants in $\mathcal{\tilde R}$ come from {numref}`rotation-matrices`:

$$
    \mathcal M=\mathcal {\tilde RH_0\tilde R}^\dagger=
    \begin{pmatrix}
        \alpha^2t_1+\beta^2t_2&\alpha\beta'(t_1-t_2)+(\alpha\alpha'-\beta\beta')t_{12}&\alpha t_{13}+\beta t_{23}\\
        \alpha\beta'(t_1-t_2)-(\alpha\alpha'-\beta\beta')t_{12}&\alpha^2t_2+\beta^2t_1&\beta't_{13}+\alpha't_{23}\\
        \alpha t_{13}-\beta t_{23}&\beta't_{13}-\alpha't_{23}&t_3
        \label{eq:rot applied}
    \end{pmatrix},
$$

Now we can express the hopping Hamiltonian as a sum of $\mathcal M$ matrices:

$$
\mathcal H_\text{hop} = \sum_{\sigma\in\{-1, 1\}}\mathcal M(1, \sigma, 0, 0)\e{i\sigma k_x}+\sum_{\sigma,\tau\in\{-1, 1\}}\mathcal M\left(-\frac12,\frac\sigma2, \frac{\tau\sqrt3}{2},\frac{\sigma\tau\sqrt3}2\right)\e{i\sigma k_x/2+i\tau k_y\sqrt3/2}
$$

````

```{list-table} Hamiltonian components with $\tilde k_x=k_x/2$ and $\tilde k_y = k_y\cdot \sqrt3/2$.
:header-rows: 1
* -  
  - Real
  - Imaginary
* - $h_1$
  - $2t_1\cos 2\tilde k_x+(t_1+3t_2)\cos\tilde k_x\cos\tilde k_y$
  - 0
* - $h_2$
  - $2t_2\cos 2\tilde k_x+(t_2+3t_1)\cos\tilde k_x\cos\tilde k_y$
  - 0
* - $h_3$
  - $2t_3\cos 2\tilde k_x+4t_3\cos\tilde k_x\cos\tilde k_y$
  - 0
* - $h_{12}$
  - $\sqrt3(t_2-t_1)\sin\tilde k_x\sin\tilde k_y$
  - $2t_{12}(\sin 2\tilde k_x-2\sin\tilde k_x\cos\tilde k_y)$
* - $h_{13}$
  - $2t_{13}(\cos2\tilde k_x-\cos\tilde k_x\cos\tilde k_y)$
  - $2\sqrt{3}t_{23}\cos\tilde k_x\sin\tilde k_y$
* - $h_{23}$
  - $-2\sqrt3t_{13}\sin\tilde k_x\sin\tilde k_y$
  - $2t_{23}(\sin 2\tilde k_x+\sin\tilde k_x\cos\tilde k_y)$
```


```{bibliography} references.bib
```

```python

```
