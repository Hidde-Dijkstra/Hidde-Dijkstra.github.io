import numpy as np
import plotly.graph_objects as go

class BandStructure:
    
    def __init__(self):
        self.fig = go.Figure()
        self.fig.update_layout(xaxis_title=r"$k$", yaxis_title= 'Energy (eV)')
    
    def plot_from_H(self, H, name, color='red', N_input=2):
        N = len(H(np.zeros(N_input)))
        bands = np.zeros((len(self.k_path), N))
        for j in range(len(self.k_path)):
            bands[j, :] = np.linalg.eigvalsh(H(self.k_path[j]))
        self.add_bands(bands, color, name)
            
    def add_bands(self, bands, color, name, k_arr=None):
        style_dict = {
                'legendgroup': name,
                'mode': 'lines',
                'line': {'color': color},
                'name': name
            }
        if np.all(k_arr == None):
            k_arr = np.arange(self.N_k)
        for i in range(0, bands.shape[1]):
            self.fig.add_trace(go.Scatter(x=k_arr, y=bands[:, i], showlegend=(i==0), **style_dict))
    
    def set_k_path(self, k_list, k_tags, N_k):
        self.N_k = N_k
        k_norms = [np.linalg.norm(k_list[i+1]-k_list[i]) for i in range(len(k_list)-1)]
        self.spacing = [int(N_k*k_norms[i]/sum(k_norms)) for i in range(len(k_norms))]
        self.k_path = []
        for i in range(len(self.spacing)):
            self.k_path += [k_list[i] + (k_list[i+1]-k_list[i])*j/self.spacing[i] for j in range(self.spacing[i])]
        self.k_path += [k_list[-1]]
        self.fig.update_xaxes(ticktext=k_tags, tickvals=[sum(self.spacing[:i]) for i in range(len(self.spacing)+1)])