import pandas as pd
import numpy as np
from webcolors import hex_to_rgb
import plotly.graph_objects as go # Import the graphical object
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot


node_label = ["A1", "A2", "B1", "B2","B3", "C1", "C2"]
node_index_dict = {'A1': 0, 'A2': 1, 'B1': 2, 'B2': 3, 'B3': 4, 'C1': 5, 'C2': 6}


source = ['A1','A1','A1','A2','A2','A2','B1','B2','B2','B3','B3']
target = ['B1','B2','B3','B1','B2','B3','C1','C1','C2','C1','C2']
values = [ 10, 5, 15, 5, 20, 45, 15, 20, 5, 30, 30 ]

source_index = [node_index_dict[x] for x in source]
target_index = [node_index_dict[x] for x in target]
# [0, 0, 0, 1, 1, 1, 2, 3, 3, 4, 4]
# [2, 3, 4, 2, 3, 4, 5, 5, 6, 5, 6]


fig = go.Figure(data=[go.Sankey(node = dict(label = node_label),
                    link = dict(source = source_index,
                                target = target_index,
                                value = values))])

# With this save the plots
plot(fig,image_filename='sankey_plot_1',image='png',image_width=1000,image_height=600)

# And shows the plot
fig.show()