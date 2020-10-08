---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.6.0
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Log Book

+++

## To-do:
* Draw lattice vectors in lattice 
* Combine svg objects into one to pass from functions so I do not need to pass the container as argument
* Investigate how to use ipywidgets to toggle between different $m$ and $n$ for bilayer lattice
* Rename files for clarity
* Figure out how to wrap function for small screen viewing (or maybe scroll bar like for tables)
* Write introduction file and update site image
* Make lattice drawing class more variable
* Include reciprocal vectors in k_path image
* Resolve reference conflicts
* Include horizontal lines in band structure
* Choose theme for plotly plots
* Write explanations for code when refactoring is done

+++

## 2-8 October

### Meeting

### 8 October

* Looking for papers on twisted bilayer WSe$_2$
    * [Deep moir´e potentials in twisted transition metal dichalcogenide bilayers](https://arxiv.org/pdf/2008.07696.pdf) about stacking types between $WSe_2$ and $MSe_2$, does not seem so useful.
    * [Atomic reconstruction and moiré patterns in transition metal
dichalcogenide van der Waals heterostructures](https://arxiv.org/pdf/1911.12282.pdf) seems to be about fractional filling of moiré supercells between MoS(e)$_2$ WS(e)$_2$. For some reason popular heterostructure.
    *

### 6-7 October

* Rewrite explanation for derivation of Hamiltonian elements.
* Change band structure code to its own class (previously part of WSe2 class).
* Figure out how to use plotly to plot functions interactively. Main issue resolved: plotting all bands for one spin to the same toggle key in the legend.

### 5 October

* Simplify symmetry explanation
* Refactoring orbital drawing code
* Include rotating orbitals for the monolayer

### 3-4 October

* Write lattice drawing class for svg lattice which also works for (rotated) bilayer lattices.  

+++

## 1 October

### Meeting 

Talked about restarting the project. Louk set out the general plan of the project and the immenent hurdles to cross:

* Correctly describe the mono-layer for both unit- and supercell
* Commentate code
* Find articles on interlayer coupling of WSe$_2$
* Consider $d$ orbital decomposition between the bilayer atoms.
