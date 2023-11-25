## Notebook 1 questions

> Why is there a drop in Gibbs free energy (a potential) after minimisation?

Energy is released into the system as a result of getting to an energy minimum.
[2nd law of thermodynamics](https://en.wikipedia.org/wiki/Second_law_of_thermodynamics)
wants to make disorder (entropy).

> How does this relate to entropy of the system?

[Chelate effect](https://en.wikipedia.org/wiki/Chelation#Chelate_effect):
entropy vs. enthalpy.

> What does `AllChem.EmbedMolecule` do?

It generates a 3D conformation of a molecule, which is stored as a `Conformer` object.
Accessible via `mol.GetConformer()`.

> What command does a 2D representation?

`AllChem.Compute2DCoords(mol)` adds 2D coordinates (via a conformer) to that molecule.

> What is the difference between a `Chem.Mol` instance and its `Chem.Conformer(s)`? And do both store cartesian atomic positions?

Chem.Mol is a molecule, Chem.Conformer is a conformation of that molecule.
The former stores atom and bond details, the latter stores coordinates.

> What is the SMILES of the compound you searched in SmallWorld? And what was the distance to a purchasable analogue?

👾👾👾

> Did you find a molecule that has a Wikipedia page but is not in EnamineREAL DB?

👾👾👾

> How many small molecules are there?

👾👾👾

> How many sites?

👾👾👾

> If you have a dimer, what do you see as a problem?

Hit binding differently in one site different due to asymmetry of the dimer.

> What data would you like to see in the above table and why?

Molecular weight, logP, number of rotatable bonds, number of hydrogen bond donors/acceptors
to get an idea of the physicochemical properties of the molecules.

Other values could be crystallographic occupancy of the hit, crystallographic blur-factors, predicted ∆G of binding, 


> What does Google say the RDKit command to do so is? (Remember than with a `pd.Series` you have the `apply` method, eg. `df.ROMol.apply(Descriptors.ExactMolWt)`

👾👾👾

> The above simply gets the molecular replacement template as the target. Is that wise?

👾👾👾

> The table has a `site_name` column. What would be a good approach to choose what sites to focus on? (Not for now, i.e. remember the adage 'a week in the lab, saves you an hour in the library').

👾👾👾

## Notebook 2 questions

> What is stored in a HETATM? (see )

👾👾👾

> What is an "apo structure"?

👾👾👾

> Why is crudely removing heteroligand atoms bad, and what could be done to fix it?

👾👾👾

> Why does the partial charge reside in an atom of a residue type not an atom type?

👾👾👾

> Why is bond order often absent in residue types/topologies?

👾👾👾

> What force dominates when atoms are too close?

👾👾👾

> What happens between 2-4 Å? (Zoom into the interactive plotly figure)

👾👾👾

## Notebook 3 questions

> What is a Monte Carlo method? Hint: it is not a method written in monégasque.