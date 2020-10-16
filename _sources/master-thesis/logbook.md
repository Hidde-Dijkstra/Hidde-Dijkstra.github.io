---
jupytext:
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

# Logbook

+++

## To-do: 
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

## 9-15 October

### Meeting

### 15 October

* Find papers on bilayer WSe$_2$:
    * [Isotope Effect in Bilayer WSe2](https://pubs.acs.org/doi/pdf/10.1021/acs.nanolett.8b04269) interlayer distance is 1.3 nm.
    * [Spin-Layer Locking Effects in Optical Orientation of Exciton Spin in
Bilayer WSe2](https://arxiv.org/pdf/1311.7087.pdf) AB stacking tight-binding useless since assumes coupling of d_$z2$ orbital to other layer is zero due to symmetry, only coupling between valence bands.

### 14 October

* Made tabbed plots for different supercells
* Wrote pseudocode to explain NN code

### 12-13 October

* Combine svg objects into groups for clearer code, apparently groups cannot be elements of groups.

### 10-11 October

* Fail to add caption to html figure, probable need to wait for update of jupyterbook

+++

## 2-8 October

### Meeting

* Find papers using people that cite relevant papers
* Consider only in plane rotations of orbitals for interlayer hopping

### 8 October

* Looking for papers on twisted bilayer WSe$_2$
    * [Deep moir´e potentials in twisted transition metal dichalcogenide bilayers](https://arxiv.org/pdf/2008.07696.pdf) about stacking types between WSe$_2$ and MSe$_2$, does not seem so useful.
    * [Atomic reconstruction and moiré patterns in transition metal
dichalcogenide van der Waals heterostructures](https://arxiv.org/pdf/1911.12282.pdf) seems to be about fractional filling of moiré supercells between MoS(e)$_2$ WS(e)$_2$. For some reason popular heterostructure.
    * [Flat bands in long moiré wavelength twisted bilayer
WSe2](https://arxiv.org/pdf/1910.13068.pdf) talks a lot about AA vs AB stacking and localization of wavefunctions but does not give information about dimensions of the lattice/coupling.
* Just do DFT ourselves?

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
