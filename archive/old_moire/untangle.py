from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
import numpy as np
import itertools

def untangle(bands, Δε, i=0):
    M, N = bands.shape
    rolls = np.abs(np.expand_dims(bands, axis=0)-np.expand_dims(bands, axis=1)) < Δε
    sorted_bands = np.zeros(bands.shape)
    sorted_bands[:, :2] = bands[:, :2]
    for n in range(2, N):
        c, cross_arr = connected_components(csr_matrix(rolls[:, :, n]), 
                                             return_labels=True)
        cross_list = [[] for _ in range(c)]
        for i, crossing in enumerate(cross_arr):
            cross_list[crossing].append(i)    
        for cross in cross_list:
            band_slice = bands[min(cross):min(cross)+len(cross), n]
            untangle_point(cross, band_slice, sorted_bands, n)
        sorted_bands = sorted_bands[np.argsort(sorted_bands[:, n])]
    return sorted_bands

def untangle_point(cross, band_slice, sorted_bands, n):
    dk2 = sorted_bands[cross, n-2] - 2*sorted_bands[cross, n-1]
    min_curvature = np.inf
    for permutation in itertools.permutations(range(len(cross))):
        curvature = np.sum(np.abs(dk2 + band_slice[list(permutation)]))
        if curvature < min_curvature:
            min_curvature = curvature
            best_permutation = list(permutation)
    sorted_bands[cross, n] = band_slice[best_permutation]
        

def untangle2(unsorted_bands, α=1):
    M, N = unsorted_bands.shape
    bands = np.zeros((M, N))
    bands[:2] = unsorted_bands[:2]
    δ_bands = bands[1] - bands[0]
    for m in range(2, M):
        δ_bands = α*(bands[m-1]-bands[m-2]) + (1-α)*δ_bands     
        bands_estimate = δ_bands + bands[m-1]
        bands = bands[:, np.argsort(bands_estimate)]
        δ_bands = δ_bands[np.argsort(bands_estimate)]
        bands[m] = unsorted_bands[m]
    return bands