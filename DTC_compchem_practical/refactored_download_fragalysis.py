"""
This is not finished! I am still working on it. I am not sure if I will be able to finish it in time.
"""


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


class FragalysisDownloader:
    fragalysis_api_url = 'fragalysis.diamond.ac.uk/api/download_structures/'
    api_data = {
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

    def __init__(self, target_name: str):
        self.target_name = target_name
        url_response: requests.Response = requests.post(f'https://{self.fragalysis_api_url}',
                                                        json={'target_name': self.target_name, **self.api_data})
        url_response.raise_for_status()
        self.file_url: str = url_response.json()['file_url']
        response: requests.Response = requests.get(f"https://{self.fragalysis_api_url}?file_url={self.file_url}",
                                                   allow_redirects=True)
        response.raise_for_status()
        self.zf = zipfile.ZipFile(io.BytesIO(response.content), "r")

    def __getitem__(self, item: str):
        for fileinfo in self.zf.infolist():
            if item in fileinfo.filename:
                return self.zf.read(fileinfo.filename).decode('utf8')

    def __iter__(self) -> Tuple[str, str]:
        for fileinfo in self.zf.infolist():
            yield fileinfo.filename, self.zf.read(fileinfo.filename).decode('utf8')

    def __len__(self):
        return len(self.zf.infolist())

    def write_all(self, directory: Optional[str] = None):
        if directory is None:
            directory = self.target_name
        if not os.path.exists(directory):
            os.makedirs(directory)
        for fileinfo in self.zf.infolist():
            if os.path.split(fileinfo.filename)[0] != '':
                os.makedirs(os.path.join(directory, os.path.split(fileinfo.filename)[0]), exist_ok=True)
            with open(os.path.join(directory, fileinfo.filename), 'w') as f:
                f.write(self.zf.read(fileinfo.filename).decode('utf8'))

    def write_key_files(self, directory: Optional[str] = None,
                        wanted_files: Sequence = ('reference.pdb', 'metadata.csv', 'combined.sdf')) -> Dict[str, str]:
        """
        Write the files of interest to a directory. set to `'.'` for not using a directory.
        """
        if directory is None:
            directory = self.target_name
        if not os.path.exists(directory):
            os.makedirs(directory)
        for wanted in wanted_files:
            self[wanted]

        for fileinfo in self.zf.infolist():

            if wanted in fileinfo.filename:
                filename = self.target_name + '_' + os.path.split(fileinfo.filename)[-1].replace(self.target_name + '_',
                                                                                                 '')
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

        for fileinfo in self.zf.infolist():
            for fileinfo in zf.infolist():
                if '.pdb' not in fileinfo.filename:
                    continue
                with open(path, 'w') as fh:
                    fh.write(zf.read(fileinfo).decode('utf8'))
                return paths
        raise ValueError('No pdb files...')

    def to_pandas(self) -> pd.DataFrame:
        # make a combined table
        # Fragalysis does not give attributes in the sdf entries. This is instead stored in metadata.csv.

        mol_df = pd.concat([PandasTools.LoadSDF(sdf_filename).set_index('ID'),
                            pd.read_csv(metadata_filename, index_col=0).set_index('crystal_name')
                            ], axis=1)

    def write_specified(self, item: str, directory: Optional[str] = None):

        output_folder: str = '.',

    paths: Dict[str, str] = {'$target_name': target_name}

