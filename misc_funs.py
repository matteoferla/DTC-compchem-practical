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
                 subImgSize=(150,150),
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
    
    
def download_fragalysis(target_name: str, 
                        output_folder: str = '.',
                        wanted_files: Sequence = ('reference.pdb', 'metadata.csv', 'combined.sdf')) -> Dict[str, str]:

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        
    api_data = {'target_name': target_name,
             'proteins': '',
             'event_info': False,
             'sigmaa_info': False,
             'diff_info': False,
             'trans_matrix_info': False,
             'NAN': False,
             'mtz_info': False,
             'cif_info': False,
             'NAN2': False,
             'map_info': False,
             'single_sdf_file': True,
             'sdf_info': False,
             'pdb_info': False,
             'bound_info': True,
             'metadata_info': True,
             'smiles_info': True,
             'static_link': False,
             'file_url': ''}

    fragalysis_api_url = 'fragalysis.diamond.ac.uk/api/download_structures/'
    response: requests.Response = requests.post(f'https://{fragalysis_api_url}', json=api_data)
    response.raise_for_status()
    file_url:str = response.json()['file_url']
    response: requests.Response = requests.get(f"https://{fragalysis_api_url}?file_url={file_url}", allow_redirects=True)
    response.raise_for_status()
    zf = zipfile.ZipFile(io.BytesIO(response.content), "r")
    paths: Dict[str, str] = {'$target_name': target_name}
    for fileinfo in zf.infolist():
        for wanted in wanted_files:
            if wanted in fileinfo.filename:
                filename = target_name+ '_' + os.path.split(fileinfo.filename)[-1].replace(target_name+'_', '')
                path = os.path.join( output_folder,  filename)
                paths[wanted] = path
                with open(path, 'w') as fh:
                    fh.write(zf.read(fileinfo).decode('utf8'))
                break
    # I am not sure why, but some lack reference.pdb
    if 'reference.pdb' not in wanted_files or 'reference.pdb' in paths:
        return paths
    path = os.path.join( output_folder,  target_name+ '_reference.pdb')
    paths['reference.pdb'] = path
    for fileinfo in zf.infolist():
        if '.pdb' not in fileinfo.filename:
            continue
        with open(path, 'w') as fh:
            fh.write(zf.read(fileinfo).decode('utf8'))
        return paths
    raise ValueError('No pdb files...')
            


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
    fig = ff.create_dendrogram(matrix, orientation='bottom', labels=mols.apply(lambda mol: mol.GetProp('_Name')).tolist(), color_threshold=1)
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
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]

    heatmap = [
        go.Heatmap(
            x = dendro_leaves,
            y = dendro_leaves,
            z = heat_data,
            colorscale = 'hot'
        )
    ]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width':1200, 'height':1200,
                             'showlegend':False, 'hovermode': 'closest',
                             })
    # Edit axes
    for domain in [[.15, 1], [0, .15], [0, .85], [.825, .975]]:
        fig.update_layout(xaxis={'domain': domain,
                                          'mirror': False,
                                          'showgrid': False,
                                          'showline': False,
                                          'zeroline': False,
                                          'ticks':""})

    return fig