from typing import *
from rdkit import Chem
import py3Dmol


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


# ================================================================
# mol grid

from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from IPython.display import display


def display_mols(mols: Sequence[Chem.Mol],
                 molsPerRow=5,
                 subImgSize=(150,150),
                 useSVG=True):
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


# ================================================================