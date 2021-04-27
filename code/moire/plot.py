import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from itertools import product
import drawSvg as draw

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
        plot = Plot(
            f, 
            cut=kwargs.get('cut'),  
            slider=kwargs.get('slider'), 
            legend=kwargs.get('legend')
        )
        plot.prepare_traces()
        return plot.plot(fig=kwargs.get('fig'))
    return wrapper

class Plot:

    default_colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    def __init__(self, f, cut=None, legend=None, slider=None, showlegend=True):
        self.showlegend = showlegend
        self.cut = cut
        self.f = f
        if legend == None:
            self.legend = dict(exists=False)
        else:
            self.legend = legend
            self.legend['exists'] = True
        if slider == None:
            self.slider = dict(exists=False)
        else:
            self.slider = slider
            self.slider['exists'] = True
            
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
        if not self.legend['exists']:
            keys.append('legend')
            self.legend['arg'] = 'legend'
            arg_list.append([None])
        for trace_id, trace in enumerated_product(*arg_list):
            f_args = {}
            f_args_id = {}
            for i, key in enumerate(keys):
                f_args[key] = trace[i]
                f_args_id[key] = trace_id[i]
            self.traces.append([f_args, f_args_id])
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
        if self.slider['exists']:
            N_traces = np.zeros(self.slider['N'], dtype=int)
        visibility_list = []
        for trace in self.traces:
            f_traces = self.f(**trace[0])
            if type(f_traces) != list:
                f_traces = [f_traces]
            for i, f_trace in enumerate(f_traces):
                visibility, visibility_list = self.check_legend_visibility(trace, visibility_list)
                f_trace.showlegend = visibility
                f_trace.legendgroup = str(trace[1][self.legend['arg']])
                if 'tags' in self.legend:
                    f_trace.name = self.legend.get('tags')[trace[1][self.legend['arg']]]
                if 'colors' in self.legend:
                    f_trace.marker['color'] = self.legend.get('colors')[trace[1][self.legend['arg']]]
                elif self.legend.get('auto_color'):
                    f_trace.marker['color'] = self.default_colors[trace[1][self.legend['arg']]%10]
            if self.cut != None:
                rows = []
                for f_trace in f_traces:
                    rows.append(1+(np.max(f_trace.y)<self.cut))
            else:
                rows = [1] * len(f_traces)
            fig.add_traces(f_traces, cols=1, rows=rows)
            if self.slider['exists']:
                N_traces[trace[1][self.slider['arg']]] += len(f_traces)

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
        i_legend = trace[1][self.legend['arg']]
        if (i_slider, i_legend) not in visibility_list:
            visibility = self.showlegend
            visibility_list.append((i_slider, i_legend))
        else:
            visibility = False
        return visibility, visibility_list

def enumerated_product(*args):
    yield from zip(product(*(range(len(x)) for x in args)), product(*args))

class BandStructure:
    
    def __init__(self, cut=None):
        self.fig = make_subplots(
            rows = 1 + (cut!=None), 
            cols = 1, 
            vertical_spacing = 0.05, 
            shared_xaxes=True,
        )
        self.cut = None
        self.fig.update_layout(xaxis_title='Momentum', yaxis_title= 'Energy (eV)')
    
    @figure
    def _plot_from_H(self, H=lambda k: 0, N_input=2, **kwargs):
        N = len(H(np.zeros(N_input)))
        bands = np.zeros((len(self.k_path), N))
        for j in range(len(self.k_path)):
            bands[j] = np.linalg.eigvalsh(H(self.k_path[j]))
        traces = []
        for band in np.transpose(bands):
            traces.append(go.Scatter(y=band))
        return traces

    def plot_from_H(self, **kwargs):
        self._plot_from_H(fig=self.fig, **kwargs)

    def add_bands(self, **kwargs):
        self._plot_bands(fig=self.fig, **kwargs)

    @figure        
    def _plot_bands(self, bands=[], **kwargs):
        traces = []
        for band in bands:
            traces.append(go.Scatter(y=band))
        return traces
    
    def set_k_path(self, k_list, k_tags, N_k):
        if type(N_k) == int:
            k_norms = [np.linalg.norm(k_list[i+1]-k_list[i]) for i in range(len(k_list)-1)]
            spacing = [int(N_k*k_norms[i]/sum(k_norms)) for i in range(len(k_norms))]
        else:
            spacing = N_k
        k_path = []
        for i in range(len(spacing)):
            k_path += [k_list[i] + (k_list[i+1]-k_list[i])*j/spacing[i] for j in range(spacing[i])]
        k_path.append(k_list[-1])
        self.k_path = k_path
        self.fig.update_xaxes(ticktext=k_tags, tickvals=[sum(spacing[:i]) for i in range(len(spacing)+1)])

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