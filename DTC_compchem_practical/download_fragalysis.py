__all__ = ['download_fragalysis']

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


def download_fragalysis(target_name: str,
                        output_folder: str = '.',
                        wanted_files: Sequence = ('reference.pdb', 'metadata.csv', 'combined.sdf')) -> Dict[str, str]:
    """
    First pass at a download function for fragalysis.
    The Fragalysis-API does it but it's rather byzantine.
    This is a bit more straightforward.
    I am refactoring this to be more object-oriented in ``refactored_download_fragalysis.py``

    :param target_name: the name of the target, e.g. 'Mpro', do note that case matters
    :param output_folder: the folder to save the files to
    :param wanted_files: the files to download, these are words in the filename.

    A special case is ``reference.pdb`` which is the reference structure for alignment.
    It is not always present, so if it is not present, the first PDB file is used instead.
    """
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
    file_url: str = response.json()['file_url']
    response: requests.Response = requests.get(f"https://{fragalysis_api_url}?file_url={file_url}",
                                               allow_redirects=True)
    response.raise_for_status()
    zf = zipfile.ZipFile(io.BytesIO(response.content), "r")
    paths: Dict[str, str] = {'$target_name': target_name}
    for fileinfo in zf.infolist():
        for wanted in wanted_files:
            if wanted in fileinfo.filename:
                filename = target_name + '_' + os.path.split(fileinfo.filename)[-1].replace(target_name + '_', '')
                path = os.path.join(output_folder, filename)
                paths[wanted] = path
                with open(path, 'w') as fh:
                    fh.write(zf.read(fileinfo).decode('utf8'))
                break
    # I am not sure why, but some lack reference.pdb
    if 'reference.pdb' not in wanted_files or 'reference.pdb' in paths:
        return paths
    path = os.path.join(output_folder, target_name + '_reference.pdb')
    paths['reference.pdb'] = path
    for fileinfo in zf.infolist():
        if '.pdb' not in fileinfo.filename:
            continue
        with open(path, 'w') as fh:
            fh.write(zf.read(fileinfo).decode('utf8'))
        return paths
    raise ValueError('No pdb files...')
