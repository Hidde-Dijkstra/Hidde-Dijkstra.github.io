from .decorators import change_directory, save, time
from . import io_processing as io
from moire.bandstructure import BandStructure, untangle
from moire.plot.figure import figure
import plotly.graph_objects as go
from sympy.physics.quantum.cg import CG
import numpy as np
import os


class PROJWFC:
    def __init__(self, DFT, projwfc_dir=""):
        self.DFT = DFT
        self.work_directory = DFT.work_directory
        self.data_directory = DFT.data_directory
        self.prefix = DFT.prefix
        self.dir = projwfc_dir

    @change_directory("work_directory")
    def write(self, ΔE=0.05):
        file = [
            f"&PROJWFC",
            f"  prefix = '{self.prefix}'",
            f"  outdir = '{self.DFT.QE.qe_dic['&CONTROL']['outdir']}'",
            f"  DeltaE = {ΔE}",
            f"/\n",
        ]
        io.write_file(f"{self.prefix}.projwfc.in", "\n".join(file))

    @change_directory("work_directory")
    @time
    def run(self, n_cores=4):
        command = (
            f"mpirun -n {n_cores} projwfc.x"
            f" < {self.prefix}.projwfc.in"
            f" > {self.prefix}.projwfc.out"
        )
        print(command)
        print("projwfc", os.system(command))

    @change_directory("data_directory")
    @save
    @change_directory("work_directory")
    def extract(self):
        f = io.read_file(f"{self.prefix}.projwfc.out")
        projwfc_data = f.split(" k =   ")
        states = io.gen_lst(projwfc_data[0], "\n     state #", seperate_state, True)
        occupations = []
        bands = []
        k_points = []
        for k_point in projwfc_data[1:]:
            band_k = []
            band_data = k_point.split("\n    |psi|^2")[:-1]
            occupation = np.zeros((len(band_data), len(states)))
            for i, band in enumerate(band_data):
                ε, ψ = band.split(" eV ====")
                band_k.append(io.scrub_str(ε.split(") = ")[1]))
                ψ = io.gen_lst(ψ, "+", lambda x: io.scrub_str(x, "*"))
                for φ in ψ:
                    occupation[i, int(φ[1]) - 1] = φ[0]
            occupations.append(occupation)
            bands.append(band_k)
            k_points.append(k_point)
        return {
            f"{self.dir}/k_points": k_points,
            f"{self.dir}/bands": np.transpose(bands),
            f"{self.dir}/states": states,
            f"{self.dir}/state_names": get_state_names(projwfc_data[0]),
            f"{self.dir}/occupations": occupations,
        }

    @change_directory("data_directory")
    def plot(
        self,
        dis_min=-np.inf,
        dis_max=np.inf,
        α=1,
        N=None,
        ΔN=None,
        print_N_wan=True,
        ε_margin=1,
        **kwargs,
    ):
        path = self.DFT.path.copy(N=N, ΔN=ΔN)
        bs = BandStructure(cut=kwargs.get("cut"), path=path)
        band_dic = dict(hoverinfo="skip", line=dict(color="black", width=0.3))
        bands = np.load(f"{self.prefix}/{self.dir}/bands.npy")
        order = untangle(bands, α=α, ignore_points=path.N_k.Σ)
        bands = bands[order]
        window = np.any(np.logical_and(dis_min < bands, bands < dis_max), axis=1)
        if print_N_wan:
            print("N wannier:", sum(window))
        ε_range = [np.min(bands[window]) - ε_margin, np.max(bands[window]) + ε_margin]
        bs.add_bands(
            bands=bands[window],
            name="bands",
            autocolor=False,
            untangle_bands=False,
            band_dic=band_dic,
            α=α,
        )
        band_dic["line"]["dash"] = "dot"
        bs.add_bands(
            bands=bands[np.logical_not(window)],
            name="bands",
            autocolor=False,
            showlegend=False,
            plot_range=ε_range,
            band_dic=band_dic,
            untangle_bands=False,
        )
        return self._plot(
            fig=bs.fig, order=order, window=window, bands=bands, path=path, **kwargs
        )

    @change_directory("data_directory")
    @figure
    def _plot(
        self,
        args_1={},
        args_2={},
        size=15,
        bands=None,
        order=None,
        window=None,
        opacity=0.7,
        **kwargs,
    ):
        path = kwargs.get("path", self.DFT.path)
        valid_states = sort_states(f"{self.prefix}/{self.dir}", args_1, args_2)
        occupations = np.load(f"{self.prefix}/{self.dir}/occupations.npy")
        densities = np.sum(
            occupations * np.expand_dims(valid_states, axis=(0, 1)), axis=2
        ).T
        return [
            go.Scattergl(
                x=[k.linear for k in path],
                y=band,
                text=density,
                opacity=opacity,
                mode="markers",
                marker=dict(size=size * density, line=dict(width=1)),
            )
            for density, band in zip(densities[order][window], bands[window])
        ]


def get_state_names(state):
    state = state.split("\n     state #")[1]
    _, q_num = state.split(", wfc")
    q_num = q_num.replace("= ", "=")
    q_num = q_num.split("(")[1].split(")")[0]
    print(["atoms"] + [q.split("=")[0] for q in q_num.split(" ")])
    return ["atoms"] + [q.split("=")[0] for q in q_num.split(" ")]


def seperate_state(state):
    atom, q_num = state.split(", wfc")
    q_num = q_num.replace("= ", "=")
    atom = int(atom.split("atom")[1].split("(")[0])
    q_num = q_num.split("(")[1].split(")")[0]
    q_num = [io.scrub_str(q) for q in q_num.split(" ")]
    return [atom] + q_num


def sort_states(prefix, args_1, args_2):
    states = np.load(f"{prefix}/states.npy")
    state_names = np.load(f"{prefix}/state_names.npy")
    spin = states.shape[1] == 4
    converted_states, key_list = convert_states(states, list(state_names), spin)
    wfc = {}
    for key in key_list:
        if key in args_1:
            wfc[key] = list_transform(args_1[key])
        elif key in args_2:
            wfc[key] = list_transform(args_2[key])
        else:
            wfc[key] = []
    return check_states(converted_states, wfc, spin)


def convert_states(states, state_names, spin):
    if spin:
        names = ["atoms", "l", "j", "m_j"]
        key_list = ["atoms", "l", "m_l", "m_s"]
    else:
        names = ["atoms", "l", "m"]
        key_list = names
    state_order = [state_names.index(name) for name in names]
    return states.T[state_order], key_list


def check_states(states, wfc, spin):
    valid_states = np.zeros(states.shape[1])
    for i, state in enumerate(states.T):
        valid_state = 0
        if spin:
            q_list = convert_q_num(*state)
        else:
            q_list = [[state, 1]]
        for q in q_list:
            valid = True
            for key, value in zip(wfc.keys(), q[0]):
                valid &= (value in wfc[key]) or (wfc[key] == [])
            valid_state += q[1] ** 2 * valid
        valid_states[i] = valid_state
    return valid_states


def convert_q_num(atoms, l, j, m_j):
    q_list = []
    for m_l in [-int(l) + n for n in range(2 * int(l) + 1)]:
        for m_s in [-0.5, 0.5]:
            q_list.append(
                [(atoms, l, m_l, m_s), float(CG(l, m_l, 0.5, m_s, j, m_j).doit())]
            )
    return q_list


def list_transform(variable):
    if type(variable) == list:
        return variable
    else:
        return [variable]
