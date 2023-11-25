import py3Dmol
def get_protein_view(pdbblock: str, resn: str, chain='A', width=800, height=400) -> py3Dmol.view:
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdbblock, 'pdb')
    view.setStyle({'chain': chain}, {'cartoon': {'color': 'turquoise'}})
    view.setStyle({'chain': {'$ne': chain}}, {})
    view.setStyle({'within':{'distance':'5', 'sel':{'resn':resn}}}, {'stick': {}})
    view.setStyle({'resn':resn}, {'stick': {'colorscheme': 'yellowCarbon'}})
    view.zoomTo({'resn':resn})
    return view

def get_mols_view(**color2mols) -> py3Dmol.view:
    view = py3Dmol.view(width=800, height=400)
    for color, mols in color2mols.items():
        for mol in mols:
            view.addModel(Chem.MolToMolBlock(mol), 'mol')
            view.setStyle({'model': -1}, {'stick': {'colorscheme': color}})
    view.zoomTo({'model': -1})
    return view