import py3Dmol
from rdkit import Chem
def get_protein_view(pdbblock: str, resn: str, chain='A', width=800, height=400) -> py3Dmol.view:
    view = py3Dmol.view(width=width, height=height)
    add_protein(view, pdbblock, resn, chain)
    view.zoomTo({'resn':resn})
    return view

def add_protein(view: py3Dmol.view, pdbblock: str, resn: str, chain='A') -> py3Dmol.view:
    view.addModel(pdbblock, 'pdb')
    view.setStyle({'chain': chain}, {'cartoon': {'color': 'turquoise'}})
    view.setStyle({'chain': {'$ne': chain}}, {})
    view.setStyle({'within':{'distance':'5', 'sel':{'resn':resn}}}, {'stick': {}})
    view.setStyle({'resn':resn}, {'stick': {'colorscheme': 'yellowCarbon'}})
    return view

def get_mols_view(**color2mols) -> py3Dmol.view:
    view = py3Dmol.view(width=800, height=400)
    add_mols(view, **color2mols)
    view.zoomTo({'model': -1})
    return view

def add_mols(view: py3Dmol.view, **color2mols) -> py3Dmol.view:
    for color, mol in color2mols.items():
        view.addModel(Chem.MolToMolBlock(mol), 'mol')
        view.setStyle({'model': -1}, {'stick': {'colorscheme': color}})
    return view