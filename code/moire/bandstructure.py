from .plot.figure import figure
from plotly.subplots import make_subplots
import numpy as np
from numpy.linalg import norm, eigvalsh
import plotly.graph_objects as go
from dataclasses import dataclass, field
from typing import List


class BandStructure:
    def __init__(self, cut=None, path=None):
        self.cut = cut
        self.path = path
        if path != None:
            self.gen_fig()

    def gen_fig(self):
        self.fig = make_subplots(
            rows=1 + (self.cut != None),
            cols=1,
            vertical_spacing=0.05,
            shared_xaxes=True,
            x_title="Momentum",
            y_title="Energy (eV)",
        )
        self.fig.update_layout(legend={"itemsizing": "constant"})
        self.fig.update_xaxes(
            ticktext=[k.tag for k in self.path.sym_points],
            tickvals=[k.linear for k in self.path.sym_points],
            constrain="domain",
        )

    def set_k_path(self, k_list, tags, N):
        self.path = Path(k_list, tags, N)
        self.gen_fig()
        return self.path

    def plot_from_H(self, name="", autocolor=True, new_plot=False, **kwargs):
        if new_plot:
            self.gen_fig()
        if "legend" not in kwargs:
            kwargs["legend"] = dict(
                args=dict(_=[None]), tags=[name], autocolor=autocolor
            )
        elif "autocolor" not in kwargs["legend"]:
            kwargs["legend"]["autocolor"] = True
        self._plot_from_H(fig=self.fig, cut=self.cut, **kwargs)

    def add_bands(self, bands, name="", autocolor=True, new_plot=False, **kwargs):
        if new_plot:
            self.gen_fig()
        if "legend" not in kwargs:
            kwargs["legend"] = dict(
                args=dict(_=[None]), tags=[name], autocolor=autocolor
            )
        self._plot_from_bands(bands=bands, fig=self.fig, cut=self.cut, **kwargs)

    @figure
    def _plot_from_H(self, H=lambda _: 0, **kwargs):
        path = kwargs.get("path", self.path)
        bands = np.transpose([eigvalsh(H(k.vector)) for k in path])
        return self._plot_bands(bands=bands, **kwargs)

    @figure
    def _plot_from_bands(self, **kwargs):
        return self._plot_bands(**kwargs)

    def _plot_bands(self, bands=[], untangle_bands=True, α=1, band_dic={}, **kwargs):
        path = kwargs.get("path", self.path)
        if untangle_bands:
            bands = bands[untangle(bands, α=α, ignore_points=path.N_k.Σ)]
        return [
            go.Scatter(x=[k.linear for k in path], y=band, **band_dic) for band in bands
        ]


@dataclass
class KPoint:
    vector: np.ndarray = np.zeros(3)
    linear: float = 0
    tag: str = ""

    def __add__(self, other):
        return KPoint(
            vector=self.vector + other.vector, linear=self.linear + other.linear
        )

    def __sub__(self, other):
        return KPoint(
            vector=self.vector - other.vector, linear=self.linear - other.linear
        )

    def __mul__(self, other):
        return KPoint(vector=self.vector * other, linear=self.linear * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return KPoint(vector=self.vector / other, linear=self.linear / other)


@dataclass
class NK:
    Σ: List[int] = field(default_factory=list)
    Δ: List[int] = field(default_factory=list)


class Path:
    def __init__(self, k_list, tags, N=None, ΔN=None):
        Δk_arr = [norm(np.subtract(k1, k2)) for k1, k2 in zip(k_list, k_list[1:])]
        self.Δk_arr = np.array(Δk_arr) / sum(Δk_arr)
        self.N_k = self.separate(N, ΔN)
        self.sym_points = [
            KPoint(vector=np.array(k), linear=k_lin, tag=tag)
            for k, k_lin, tag in zip(
                k_list, np.cumsum(np.pad(self.Δk_arr, (1, 0))), tags
            )
        ]
        self.path = [
            k1 + (k2 - k1) * j / Δn
            for k1, k2, Δn in zip(self.sym_points, self.sym_points[1:], self.N_k.Δ)
            for j in range(Δn)
        ]
        self.path.append(self.sym_points[-1])

    def __getitem__(self, i=0):
        return self.path[i]

    def __len__(self):
        return len(self.path)

    def copy(self, N=None, ΔN=None):
        if np.all(N == None) and np.all(ΔN == None):
            return self
        else:
            return Path(
                k_list=[k.vector for k in self.sym_points],
                tags=[k.tag for k in self.sym_points],
                N=N,
                ΔN=ΔN,
            )

    def separate(self, N, ΔN):
        if type(N) == int:
            N -= 1
            δn = np.array([int(N * Δk) for Δk in self.Δk_arr], dtype=int)
            δN = N - np.sum(δn)
            if δN != 0:
                rest = self.Δk_arr - δn / N
                δn[rest.argsort()[-δN:]] += 1
            ΔN = δn
        elif np.all(N != None):
            return NK(Σ=N, Δ=np.array([n2 - n1 for n1, n2 in zip(N, N[1:])]))
        return NK(Σ=np.cumsum(np.pad(ΔN, (1, 0))), Δ=ΔN)


def untangle(unsorted_bands, α=1, **kwargs):
    bands = np.sort(unsorted_bands, axis=0)
    N, M = bands.shape
    ignore_points = kwargs.get("ignore_points", [0, M - 1])
    band_slices = np.zeros((M, N))
    index = np.zeros((M, N), dtype=int)
    for i, point in enumerate(ignore_points[:-1]):
        index[point] = range(N)
        index[point + 1] = range(N)
        band_slices[point : point + 2] = bands[:, point : point + 2].T
        δ_bands = band_slices[point + 1] - band_slices[point]
        for m in range(2 + point, ignore_points[i + 1] + 1):
            δ_bands = α * (band_slices[m - 1] - band_slices[m - 2]) + (1 - α) * δ_bands
            bands_estimate = δ_bands + band_slices[m - 1]
            band_slices = band_slices[:, np.argsort(bands_estimate)]
            δ_bands = δ_bands[np.argsort(bands_estimate)]
            index = index[:, np.argsort(bands_estimate)]
            band_slices[m] = bands[:, m]
            index[m] = range(N)
    return index.T, np.arange(M)[:, None].T
