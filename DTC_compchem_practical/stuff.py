from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools, rdShapeHelpers, rdmolops, Descriptors
from IPython.display import display
from typing import Sequence, Optional, List, Dict, Tuple
import requests, json, zipfile, io, os

import plotly.graph_objects as go
import plotly.figure_factory as ff
import pandas as pd

import numpy as np
from scipy.spatial.distance import pdist, squareform


def display_mols(mols: Sequence[Chem.Mol],
                 molsPerRow=5,
                 subImgSize=(150, 150),
                 useSVG=True) -> None:
    """
    Rudimentary wrapper for calling ``display(Draw.MolsToGridImage``
    """
    flattos = [AllChem.RemoveHs(mol) for mol in mols if isinstance(mol, Chem.Mol)]
    for mol in flattos:
        AllChem.Compute2DCoords(mol)
    display(Draw.MolsToGridImage(flattos,
                                 legends=[mol.GetProp('_Name') if mol.HasProp('_Name') else '-' for mol in mols],
                                 subImgSize=subImgSize, useSVG=useSVG,
                                 molsPerRow=molsPerRow))


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


def distance_heatmap(df) -> go.Figure:
    # all of this stolen from Rachael, who stole it straight from plotly: https://plotly.com/python/dendrogram/

    # distance
    mols = df.ROMol.apply(rdmolops.AddHs)
    matrix: np.array = calc_distances(mols)

    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(matrix, orientation='bottom',
                               labels=mols.apply(lambda mol: mol.GetProp('_Name')).tolist(), color_threshold=1)
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
    for domain in [[.15, 1], [0, .15], [0, .85], [.825, .975]]:
        fig.update_layout(xaxis={'domain': domain,
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'ticks': ""})

    return fig