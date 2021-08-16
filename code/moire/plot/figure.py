import numpy as np
from plotly.subplots import make_subplots
from itertools import product
import plotly


def figure(func):
    def wrapper(self=None, **kwargs):
        def f(**kwargs2):
            return func(self=self, **kwargs, **kwargs2)

        plot = Plot(f, **kwargs)
        plot.prepare_traces()
        return plot.plot(fig=kwargs.get("fig"))

    return wrapper


class Plot:

    default_colors = [
        "#636EFA",
        "#EF553B",
        "#00CC96",
        "#AB63FA",
        "#FFA15A",
        "#19D3F3",
        "#FF6692",
        "#B6E880",
        "#FF97FF",
        "#FECB52",
    ]
    default_markers = ["circle", "square", "diamond", "cross", "x", "star"]
    i = 0

    def __init__(
        self,
        f,
        showlegend=True,
        autolegend=True,
        cut=None,
        plot_range=[-np.inf, np.inf],
        reset_color=False,
        **kwargs
    ):
        if kwargs.get("legend") == None:
            self.legend = dict(exists=False)
        else:
            self.legend = kwargs.get("legend")
            self.legend["exists"] = True
        if kwargs.get("slider") == None:
            self.slider = dict(exists=False)
        else:
            self.slider = kwargs.get("slider")
            self.slider["exists"] = True
        self.plot_range = np.array(plot_range)
        self.autolegend = autolegend
        self.showlegend = showlegend
        self.cut = cut
        self.f = f
        if reset_color:
            Plot.i = 0

    def prepare_traces(self):
        self.traces = []
        keys = []
        arg_list = []
        for dic in [self.slider, self.legend]:
            if dic["exists"]:
                key, args = next(iter(dic.get("args").items()))
                keys.append(key)
                arg_list.append(args)
                dic["arg"] = key
                dic["N"] = len(args)
        for trace_id, trace in enumerated_product(*arg_list):
            f_args = {}
            f_args_id = {}
            for i, key in enumerate(keys):
                f_args[key] = trace[i]
                f_args_id[key] = trace_id[i]
            self.traces.append([f_args, f_args_id])
        if arg_list == []:
            self.traces.append([{}, {}])
        if self.slider["exists"]:
            self.traces.sort(key=lambda x: x[1].get(self.slider["arg"]))

    def plot(self, fig=None):
        if fig == None:
            fig = make_subplots(
                rows=1 + (self.cut != None),
                cols=1,
                vertical_spacing=0.05,
                shared_xaxes=True,
            )
        N_0 = len(fig.data)
        Δi = {}
        if self.slider["exists"]:
            N_traces = np.zeros(self.slider["N"], dtype=int)
        visibility_list = []
        for trace in self.traces:
            f_traces = self.f(**trace[0])
            if type(f_traces) != list:
                if type(f_traces) == plotly.graph_objs._figure.Figure:
                    fig.update_layout(f_traces.layout)
                    f_traces = list(f_traces.data)
                else:
                    raise Exception("invalid input")
            rows = []
            valid_traces = []
            for f_trace in list_transform(f_traces):
                range_bool, row = self.in_range(f_trace.y)
                if range_bool:
                    valid_traces.append(f_trace)
                    rows.append(row)
            for i, f_trace in enumerate(valid_traces):
                visibility, visibility_list = self.check_legend_visibility(
                    trace, visibility_list
                )
                if self.autolegend:
                    f_trace.showlegend = visibility
                if self.legend["exists"]:
                    f_trace.legendgroup = str(
                        trace[1][self.legend["arg"]]
                    )  # +str(self.show_legend)
                    if "tags" in self.legend:
                        f_trace.name = self.legend.get("tags")[
                            trace[1][self.legend["arg"]]
                        ]
                        f_trace.legendgroup += f_trace.name
                    if "colors" in self.legend:
                        f_trace.marker["color"] = self.legend.get("colors")[
                            trace[1][self.legend["arg"]]
                        ]
                    elif self.legend.get("autocolor"):
                        i = trace[1][self.legend["arg"]]
                        f_trace.marker["color"] = self.default_colors[(Plot.i + i) % 10]
                        Δi[str(i)] = None
                    if self.legend.get("autosymbol"):
                        f_trace.marker["symbol"] = (
                            self.default_markers[trace[1][self.legend["arg"]] % 6]
                            + "-open"
                        )
            fig.add_traces(valid_traces, cols=1, rows=rows)
            if self.slider["exists"]:
                N_traces[trace[1][self.slider["arg"]]] += len(valid_traces)
        Plot.i += len(Δi)
        if self.slider["exists"]:
            steps = []
            for i in range(self.slider["N"]):
                if "tags" in self.slider:
                    label = self.slider["tags"][i]
                else:
                    label = str(i)
                if "prefix" in self.slider:
                    currentvalue = {"prefix": self.slider["prefix"] + ": "}
                else:
                    currentvalue = {}
                step = dict(
                    method="update",
                    args=[
                        dict(visible=[True] * N_0 + [False] * np.sum(N_traces))
                    ],  # layout attribute
                    label=label,
                )
                for j in range(N_traces[i]):
                    step["args"][0]["visible"][N_0 + np.sum(N_traces[:i]) + j] = True
                steps.append(step)
            if "active" in self.slider:
                active = self.slider["active"]
            else:
                active = 0
            for i in range(self.slider["N"]):
                for j in range(N_traces[i]):
                    fig.data[N_0 + np.sum(N_traces[:i]) + j].visible = i == active
            fig.update_layout(
                sliders=[dict(steps=steps, active=active, currentvalue=currentvalue)]
            )
        return fig

    def check_legend_visibility(self, trace, visibility_list):
        if self.slider["exists"]:
            i_slider = trace[1][self.slider["arg"]]
        else:
            i_slider = 0
        if self.legend["exists"]:
            i_legend = trace[1][self.legend["arg"]]
        else:
            i_legend = 0
        if (i_slider, i_legend) not in visibility_list:
            visibility = self.showlegend
            visibility_list.append((i_slider, i_legend))
        else:
            visibility = False
        return visibility, visibility_list

    def in_range(self, arr):
        if self.cut == None:
            return np.all((self.plot_range[0] < arr) & (self.plot_range[1] > arr)), 1
        else:
            if np.all((self.plot_range[0] < arr) & (self.cut > arr)):
                return True, 2
            elif np.all((self.plot_range[1] > arr) & (self.cut < arr)):
                return True, 1
            else:
                return False, 1


def enumerated_product(*args):
    yield from zip(product(*(range(len(x)) for x in args)), product(*args))


def list_transform(variable):
    if type(variable) == list:
        return variable
    else:
        return [variable]
