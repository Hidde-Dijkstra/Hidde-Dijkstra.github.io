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

# Bloch Theorem and Tight-Binding Model

+++

We study the behavior of electrons in periodic lattices of atoms. Every atom is identical and as such we have translational invariance for any $\mathbf T=\sum n_i\mathbf a_i$, where the $a_i$ are the lattice vectors of the lattice. We model the lattice with a periodic potential $V(\mathbf r)$ that obeys the translational symmetry: $V(\mathbf r)=V(\mathbf r+\mathbf T)$. Now the stationary Schrödinger equation describes the behavior of an electron in this potential:
```{math}
H=-\frac{1}{2m}\nabla_{\mathbf r}^2+V(\mathbf r),
```
here, and in the future, we take $\hbar =1$.

```{code-cell} ipython3

```