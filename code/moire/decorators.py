import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from itertools import product

def my_plot(func):
    decorators = [layout, slider, legend, cut]
    for decorator in reversed(decorators):
        func = decorator(func)
    return func


def layout(func):
    def wrapper(**kwargs):
        return func(**kwargs)
    return wrapper

def slider(func):
    def wrapper(**kwargs):
        if 'cut' in kwargs:
            N_subplots = 2
            fig = make_subplots(
                rows=2, cols=1, 
                vertical_spacing = 0.05, 
                shared_xaxes=True,
            )
        else:
            N_subplots = 1
            fig = go.Figure()
        if 'slider' in kwargs:
            slider = kwargs['slider']
            key, args = list(slider['args'].items())[0]
            N_slider = len(args)
            for i in range(len(args)):
                traces = func(**{key: args[i]}, **kwargs)
                N_traces = 1
                if type(traces) == list:
                    N_legends = len(traces)
                    for trace in traces:
                        for i, subplots in enumearte(trace):
                            if type(subplots) == list:
                                for subplot in subplots:
                                    if type(subplot) == list:
                                        for band in subplot:
                                            fig.add_trace(band, row=i)
                            else:
                                pass
                else:
                    fig.add_trace(traces)
            steps = []
            N_traces = N_legends * N_subplots * N_trace
            for i in range(N_slider):
                step = dict(
                    method="update",
                    args=[{"visible": [False] * len(fig.data)}],  # layout attribute
                    label=str(i)
                )
                for j in range(N_traces):
                    step['args'][0]['visible'][N_traces*i+j] = True
                steps.append(step)
            if 'active' in slider:
                active = slider['active']
            else:
                active = 0
            for i in range(N_traces):
                fig.data[active+i].visible = True
            fig.update_layout(
                sliders = [dict(steps=steps, active=active)]
            )
            return fig
        else:
            return fig.add_trace(func(**kwargs))
    return wrapper


def legend(func):
    def wrapper(**kwargs):
        if 'legend' in kwargs:
            legend = kwargs['legend']
            key, args = list(legend['args'].items())[0]
            traces = []
            for i in range(len(args)):
                traces.append(func(**{key: args[i]}, **kwargs))
            return traces
        else:
            return func(**kwargs)
    return wrapper

def cut(func):
    def wrapper(**kwargs):
        if 'cut' in kwargs:
            output = func(**kwargs)
            N_cut = np.max(np.sum(output>cut, axis=0))
            return [output[:N_cut, :], output[N_cut:, :]]
        else:
            return [func(kwargs)]
    return wrapper

class Plot:
    
    def __init__(self, f):
        slider = False
        legend = False
        cut = False
        self.f = f
        self.dic_list = []
        
    def set_legend(self, legend_dic):
        if legend_dic != None:
            self.dict_list.append(legend_dic)
            self.legend = True
            
    def prepare_traces(self):
        traces = []
        keys = []
        arg_list = []
        for dic in self.dic_list:
            key, args = next(iter(dic.items()))
            keys.append(key)
            arg_list.append(args)
        for trace_id, trace in enumerated_product(arg_list):
            f_args = {}
            f_args_id = {}
            for i, key in enumerate(keys):
                f_args[key] = trace[i]
                f_args_id[key] = trace_id[i]
            traces.append([f_args, f_args_id])
        if self.slider:
            traces.sort(key=lambda x: x[1].get('slider'))

for x in itertools.product(*arg_list):
    kwargs = {}
    for i, key in enumerate(keys):
        kwargs[key] = x[i]
    print(kwargs)
    
    
def enumerated_product(*args):
    yield from zip(product(*(range(len(x)) for x in args)), product(*args))