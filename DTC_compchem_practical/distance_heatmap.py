""""
In this submodule there are two main functions.

``calc_distance_heatmap`` is a refactoring of ``original``,
but is broken...

"""

from typing import Dict, Tuple, List
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers, rdmolops, Descriptors
from scipy.spatial.distance import pdist, squareform

def calc_distances(mols: pd.Series) -> np.array:
    # calculate distance matrix
    dejavu: Dict[Tuple, float] = {}

    def calc_d(mol1: Chem.Mol, mol2: Chem.Mol):
        name_pair = tuple(sorted([mol1.GetProp('_Name'), mol2.GetProp('_Name')]))
        if name_pair in dejavu:
            return dejavu[name_pair]
        light, heavy = (mol2, mol1) if Descriptors.ExactMolWt(mol2) < Descriptors.ExactMolWt(mol1) else (mol1, mol2)
        d = rdShapeHelpers.ShapeProtrudeDist(light, heavy, allowReordering=True, ignoreHs=True)
        dejavu[name_pair] = d
        return d

    return np.array([mols.apply(lambda mol2: calc_d(mol, mol2)).tolist() for mol in mols])


def calc_distance_heatmap(mols: List[Chem.Mol]) -> go.Figure:
    # all of this stolen from Rachael, who stole it straight from plotly: https://plotly.com/python/dendrogram/

    # distance
    mol_series = pd.Series({mol.GetProp('_Name'): mol for mol in mols})
    calc_distance_heatmap()
    hmols = mol_series.apply(rdmolops.AddHs)
    matrix: np.array = calc_distances(hmols)

    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(matrix, orientation='bottom',
                               labels=hmols.index.tolist(), color_threshold=1)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(matrix, orientation='right', color_threshold=1)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(matrix)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = [
        go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=heat_data,
            colorscale='hot'
        )
    ]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width': 1200, 'height': 1200,
                       'showlegend': False, 'hovermode': 'closest',
                       })
    # Edit axes
    for axis, domain in dict(xaxis=[.15, 1],
                             xaxis2=[0, .15],
                             yaxis=[0, .85],
                             yaxis2=[.825, .975]):
        fig.update_layout(**{axis: {'domain': domain,
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'ticks': ""}})

    return fig

# ------------------------------ original code ------------------------------

import plotly.graph_objects as go

import numpy as np
from scipy.spatial.distance import pdist, squareform

def original():
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(full_matrix, orientation='bottom', labels=list(df['ID']), color_threshold=1)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(full_matrix, orientation='right', color_threshold=1)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(full_matrix)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]

    heatmap = [
        go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=heat_data,
            colorscale='hot'
        )
    ]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width': 1200, 'height': 1200,
                       'showlegend': False, 'hovermode': 'closest',
                       })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [.15, 1],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'ticks': ""})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .15],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ""})

    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'showticklabels': False,
                             'ticks': ""
                             })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain': [.825, .975],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ""})

    # Plot!
    fig.show()