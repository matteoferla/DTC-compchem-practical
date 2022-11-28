__all__ = ['make_3Dview', 'pose_3Dview', 'stylize']

"""
This submodule is to use py3Dmol to display molecules in Jupyter notebooks.
But I prefer NGLView, which does not work well in Colab...
"""

# import py3Dmol
# view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
# view.addModel(Chem.MolToMolBlock(mol2), "mol")
# display(view)

from typing import *
from rdkit import Chem
import py3Dmol
import pyrosetta_help as ph


def stylize(representation: str, color: str) -> Dict[str, Dict[str, Dict[str, str]]]:
    if 'carbon' in color.lower():
        return dict(style={representation: {'colorscheme': color}})
    else:
        return dict(style={representation: {'color': color}})


def make_3Dview(template_pdbblock, colormols: Dict[str, List[Chem.Mol]]) -> py3Dmol.view:
    """
    colormols is a diction of color/colorscheme to list of mols.
    """
    view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
    view.addModel(template_pdbblock, "pdb", stylize('cartoon', 'gainsboro'))
    view.setStyle(dict(hetflag=True), stylize('stick', 'whiteCarbon'))
    for color, mols in colormols.items():
        for mol in mols:
            view.addModel(Chem.MolToMolBlock(mol), "mol", stylize('stick', color))
    view.zoomTo(dict(hetflag=True))
    return view

def pose_3Dview(pose) -> py3Dmol.view:
    view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
    view.addModel(ph.get_pdbstr(pose), "pdb", stylize('cartoon', 'gainsboro'))
    view.setStyle(dict(hetflag=True), stylize('stick', 'whiteCarbon'))
    view.zoomTo(dict(hetflag=True))
    return view