from typing import Sequence
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display


def display_mols(mols: Sequence[Chem.Mol],
                 molsPerRow=5,
                 subImgSize=(150, 150),
                 useSVG=True) -> None:
    """
    Rudimentary wrapper for calling ``display(Draw.MolsToGridImage(...))``
    """
    flattos = [AllChem.RemoveHs(mol) for mol in mols if isinstance(mol, Chem.Mol)]
    for mol in flattos:
        AllChem.Compute2DCoords(mol)
    display(Draw.MolsToGridImage(flattos,
                                 legends=[mol.GetProp('_Name') if mol.HasProp('_Name') else '-' for mol in mols],
                                 subImgSize=subImgSize, useSVG=useSVG,
                                 molsPerRow=molsPerRow))