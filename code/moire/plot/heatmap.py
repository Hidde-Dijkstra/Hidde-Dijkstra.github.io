from scipy.interpolate import griddata
from itertools import product
import numpy as np
import plotly.graph_objs as go
from .figure import figure

a_1 = np.array([1, 0, 0])
a_2 = np.array([-1/2, np.sqrt(3)/2, 0])

def gen_heatmap(a=1, **kwargs):
    fig = _gen_heatmap(a=a, showlegend=False, **kwargs)
    fig.update_layout(xaxis=dict(title='Å', constrain='domain'),
                      yaxis=dict(scaleanchor="x", title='Å',
                                 constrain='domain', scaleratio=1))
    return fig

@figure
def _gen_heatmap(positions=[], a=1, ΔR=0.2, N=100, plot_text=False, **kwargs):
    text = go.Scatter(x=np.array([0, 1/2, 0])*a,
                      y=np.array([0, np.sqrt(3)/6, np.sqrt(3)/3])*a,
                      mode="text",
                      textfont=dict(size=20, color='white'),
                      text=["MM", "XX", "XM"])
    x, y, z = np.transpose([r+a*a_1*i+a*a_2*j for i, j, r 
                            in product(range(-2, 3), range(-2, 3), positions)])
    xi, yi = np.meshgrid(np.linspace(-ΔR+a_2[0], 1+ΔR, N), 
                         np.linspace(-ΔR, ΔR+a_2[1], N))
    heatmap = go.Heatmap(x=a*np.linspace(-ΔR+a_2[0], 1+ΔR, N),
                         y=a*np.linspace(-ΔR, ΔR+a_2[1], N),
                         z=griddata(np.transpose([x, y]), z, (a*xi, a*yi)))
    if plot_text:
        return [heatmap, text]
    else:
        return [heatmap]