import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from itertools import product
import drawSvg as draw
import plotly

def arrow(start, end, stroke_width=0.1, stroke='black', **kwargs):
    start, end = np.array(start), np.array(end)
    Δx = 3
    my_arrow = draw.Marker(-1+Δx/4, -0.5, Δx/4, 0.5, scale=4, orient='auto')
    my_arrow.append(draw.Lines(-1+Δx/4, -0.5, -1+Δx/4, 0.5, Δx/4, 0, close=True, fill=stroke))
    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none',
              marker_end=my_arrow, **kwargs)
    t = 1 - stroke_width*Δx/np.linalg.norm(end-start)
    return p.M(*start).L(*(t*(end-start)+start))

def figure(func):
    def wrapper(self=None, **kwargs):
        def f(**kwargs2):
            return func(self=self, **kwargs, **kwargs2)
        plot = Plot(f, **kwargs)
        plot.prepare_traces()
        return plot.plot(fig=kwargs.get('fig'))
    return wrapper

class Plot:

    default_colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    default_markers = ['circle', 'square', 'diamond', 'cross', 'x', 'star']
    i = 0

    def __init__(self, f, showlegend=True, autolegend=True, cut=None, 
                 plot_range=[-np.inf, np.inf], reset_color=False, **kwargs):
        if kwargs.get('legend') == None:
            self.legend = dict(exists=False)
        else:
            self.legend = kwargs.get('legend')
            self.legend['exists'] = True
        if kwargs.get('slider') == None:
            self.slider = dict(exists=False)
        else:
            self.slider = kwargs.get('slider')
            self.slider['exists'] = True
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
            if dic['exists']:
                key, args = next(iter(dic.get('args').items()))
                keys.append(key)
                arg_list.append(args)
                dic['arg'] = key
                dic['N'] = len(args)
        for trace_id, trace in enumerated_product(*arg_list):
            f_args = {}
            f_args_id = {}
            for i, key in enumerate(keys):
                f_args[key] = trace[i]
                f_args_id[key] = trace_id[i]
            self.traces.append([f_args, f_args_id])
        if arg_list == []:
            self.traces.append([{}, {}])
        if self.slider['exists']:
            self.traces.sort(key=lambda x: x[1].get(self.slider['arg']))

    def plot(self, fig=None):
        if fig == None:
            fig = make_subplots(
                rows = 1 + (self.cut!=None), 
                cols = 1, 
                vertical_spacing = 0.05, 
                shared_xaxes=True,
            )
        N_0 = len(fig.data)
        Δi = {}
        if self.slider['exists']:
            N_traces = np.zeros(self.slider['N'], dtype=int)
        visibility_list = []
        for trace in self.traces:
            f_traces = self.f(**trace[0])
            if type(f_traces) != list:
                if type(f_traces) == plotly.graph_objs._figure.Figure:
                    fig.update_layout(f_traces.layout)
                    f_traces = list(f_traces.data)
                else:
                    raise Exception('invalid input')
            rows = []
            valid_traces = []
            for f_trace in list_transform(f_traces):
                range_bool, row = self.in_range(f_trace.y)
                if range_bool:
                    valid_traces.append(f_trace)
                    rows.append(row)
            for i, f_trace in enumerate(valid_traces):
                visibility, visibility_list = self.check_legend_visibility(trace, visibility_list)
                if self.autolegend:
                    f_trace.showlegend = visibility
                if self.legend['exists']:
                    f_trace.legendgroup = str(trace[1][self.legend['arg']])#+str(self.show_legend)
                    if 'tags' in self.legend:
                        f_trace.name = self.legend.get('tags')[trace[1][self.legend['arg']]]
                        f_trace.legendgroup += f_trace.name
                    if 'colors' in self.legend:
                        f_trace.marker['color'] = self.legend.get('colors')[trace[1][self.legend['arg']]]
                    elif self.legend.get('autocolor'):
                        i = trace[1][self.legend['arg']]
                        f_trace.marker['color'] = self.default_colors[(Plot.i+i)%10]
                        Δi[str(i)] = None
                    if self.legend.get('autosymbol'):
                        f_trace.marker['symbol'] = self.default_markers[trace[1][self.legend['arg']]%6]+'-open'
            fig.add_traces(valid_traces, cols=1, rows=rows)
            if self.slider['exists']:
                N_traces[trace[1][self.slider['arg']]] += len(valid_traces)
        Plot.i += len(Δi)
        if self.slider['exists']:
            steps = []
            for i in range(self.slider['N']):
                if 'tags' in self.slider:
                    label = self.slider['tags'][i]
                else:
                    label = str(i)
                if 'prefix' in self.slider:
                    currentvalue = {'prefix': self.slider['prefix']+': '}
                else:
                    currentvalue = {}
                step = dict(
                    method="update",
                    args=[dict(visible=[True]*N_0+[False]*np.sum(N_traces))],  # layout attribute
                    label=label
                )
                for j in range(N_traces[i]):
                    step['args'][0]['visible'][N_0+np.sum(N_traces[:i])+j] = True
                steps.append(step)
            if 'active' in self.slider:
                active = self.slider['active']
            else:
                active = 0
            for i in range(self.slider['N']):
                for j in range(N_traces[i]):
                    fig.data[N_0+np.sum(N_traces[:i])+j].visible = (i==active)
            fig.update_layout(
                sliders = [dict(steps=steps, active=active, currentvalue=currentvalue)]
            )
        return fig
    
    def check_legend_visibility(self, trace, visibility_list):
        if self.slider['exists']:
            i_slider = trace[1][self.slider['arg']]
        else:
            i_slider = 0
        if self.legend['exists']:
            i_legend = trace[1][self.legend['arg']]
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
            return np.all((self.plot_range[0]<arr)&(self.plot_range[1]>arr)), 1
        else:
            if np.all((self.plot_range[0]<arr)&(self.cut>arr)):
                return True, 2
            elif np.all((self.plot_range[1]>arr) & (self.cut<arr)):
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

class BandStructure:
    
    def __init__(self, cut=None):
        self.fig = make_subplots(
            rows = 1 + (cut!=None), 
            cols = 1, 
            vertical_spacing = 0.05, 
            shared_xaxes=True,
            x_title='Momentum', 
            y_title= 'Energy (eV)',
        )
        self.fig.update_layout(legend= {'itemsizing': 'constant'})
        self.cut = cut
        dic = dict(
            legendgroup='none',
            hoverinfo='skip',
            showlegend=False, 
            line=dict(color='blue'))
    
    def plot_from_H(self, name='', showlegend=True, autocolor=True, **kwargs):
        legend = dict(args=dict(_=[None]), tags=[name], autocolor=autocolor)
        self._plot_from_H(fig=self.fig, cut=self.cut, legend=legend, showlegend=showlegend, **kwargs)

    def add_bands(self, bands, name='', showlegend=True, autocolor=True, k_arr=None, **kwargs):
        legend = dict(args=dict(_=[None]), tags=[name], autocolor=autocolor)
        self._plot_bands(bands=bands, fig=self.fig, cut=self.cut, legend=legend, showlegend=showlegend,
            k_arr=k_arr, **kwargs)
    @figure
    def _plot_from_H(self, H=lambda k: 0, N_input=2, k_arr=None, **kwargs):
        N = len(H(np.zeros(N_input)))
        bands = np.zeros((len(self.k_path), N))
        if 'band_dic' in kwargs:
            band_dic = kwargs['band_dic']
        else:
            band_dic = {}
        for j in range(len(self.k_path)):
            bands[j] = np.linalg.eigvalsh(H(self.k_path[j]))
        traces = []
        for band in np.transpose(bands):
            traces.append(go.Scatter(x=self.k_spacing, y=band, **band_dic))
        return traces

    @figure        
    def _plot_bands(self, bands=[], k_arr=None, **kwargs):
        if 'band_dic' in kwargs:
            band_dic = kwargs['band_dic']
        else:
            band_dic = {}
        traces = []
        for band in bands:
            if np.all(k_arr == None):
                traces.append(go.Scatter(x=self.k_spacing, y=band, **band_dic))
            else:
                traces.append(go.Scatter(x=k_arr, y=band, **band_dic))
        return traces
    
    def set_k_path(self, k_list, k_tags, N_k):
        k_list = [np.array(k) for k in k_list]
        k_norms = [np.linalg.norm(k_list[i+1]-k_list[i]) for i in range(len(k_list)-1)]
        if type(N_k) == int:
            n = [int(N_k*k_norms[i]/sum(k_norms)) for i in range(len(k_norms))] 
            ΔN_k = N_k - sum(n)
            if ΔN_k != 0:
                rest = np.array(k_norms) - sum(k_norms) * np.array(n)/N_k
                np.array(n, dtype=int)[rest.argsort()[-ΔN_k:]] += 1
            self.spacing = n
        else:
            self.spacing  = N_k
        k_path = []
        k_spacing = []
        for i, Δn in enumerate(self.spacing):
            k_path += [k_list[i] + (k_list[i+1]-k_list[i])*j/Δn for j in range(Δn)]
            k_spacing += [sum(k_norms[:i]) + (k_norms[i])*j/Δn for j in range(Δn)]
        k_path.append(k_list[-1])
        k_spacing.append(sum(k_norms))
        self.k_spacing = np.array(k_spacing) / sum(k_norms)
        self.k_path = k_path
        tickvals = [self.k_spacing[sum(self.spacing[:i])] for i in range(len(self.spacing)+1)]
        self.fig.update_xaxes(ticktext=k_tags, tickvals=tickvals)

class Orbital:  
    
    
    def lobe(self, color, rotate=0, translate=(0, 0), stroke="black", **kwargs):
        gradient = draw.RadialGradient(0, 1, 0.5)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(np.sqrt(3), color, 0.7)
        transform = "translate(" + " ".join([str(i) for i in translate]) + ")\nrotate(" + str(rotate) + " 0 0)"
        my_path = "M 0,0 C " + str(-np.sqrt(3)) + ",-2 " + str(np.sqrt(3)) +",-2 0,0 z"
        return draw.Path(d=my_path, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs)
    
    def circle(self, color, rotate=0, translate=(0, 0), stroke="black", ellipse=False, **kwargs):
        gradient = draw.RadialGradient(0, 0, 0.5)
        gradient.addStop(0, 'white', 0.7)
        gradient.addStop(np.sqrt(3), color, 0.7)
        transform = "rotate(" + str(rotate) + " 0 0)\ntranslate(" + " ".join([str(i) for i in translate]) + ")"
        if ellipse:
            clip = draw.ClipPath()
            clip.append(draw.Ellipse(0, 0, 0.5, 0.125, transform=transform))
            return draw.Ellipse(0, 0, 1, 0.25, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs) 
        else:
            return draw.Circle(0, 0, 0.5, stroke=stroke, stroke_width=0.01, fill=gradient, transform=transform, **kwargs)
    
    def d_xy(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
             **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=85+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=95+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=275+rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=265+rotate, translate=translate))
        return group
    
    def d_z2(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
             **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=180+rotate, translate=translate))
        group.append(self.circle(pos_color, ellipse=True, rotate=rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=rotate, translate=translate))
        return group
    
    def d_x2y2(self, translate=(0, 0), rotate=0, neg_color="dodgerblue", pos_color="red",
               **kwargs):
        group = draw.Group(**kwargs)
        group.append(self.lobe(neg_color, rotate=180+rotate, translate=translate))
        group.append(self.lobe(neg_color, rotate=rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=90+rotate, translate=translate))
        group.append(self.lobe(pos_color, rotate=270+rotate, translate=translate))
        return group