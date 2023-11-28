# Notebook 1

### Questions in first cell

> This is an enzyme, which binds to a nucleotide modification on protein. 
> Can you guess where its native substrate binds? 

In the thermodynamic sink where the bulk of the chemical matter binds

> Which fragment site would be good for a competitive inhibitor?

Active site

> Do all fragments bind to the protein with the same count of contacts per atom? If not, what is ligand efficiency?

No. Ligand efficiency is the ratio of the binding affinity / Gibbs to the number of heavy atoms.

>  Is the protein rigid? I.e. lock & key Vs. induced fit.

Flexible

> Towards the end a loop "everts", reflecting a native product bound state vs. a native substrate bound state. Would repacking sidechains model this change?

Backbone movement

> Would you reckon, the apo structure choice affect docking results?

Template choice is important.

> Some compounds have labels like `ZINC922`. What is Zinc?

A database of compounds (https://zinc20.docking.org/substances/home/) by John Irwin.

> The theoretical merger has a ester bond, while the analogue in make-on-demand space has an amide. Can you guess why?

Synthetic accessibility

> What makes a lead, "drug-like"? (cf. [Lipinski rule_of 5](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five) )
> And is a carboxyl group good? What is a bioisostere?

Crosses membranes and does not get metabolised. Carboxyl group is not good as it is charged at physiological pH,
so cannot cross membranes. 
A bioisostere is a molecule that has similar physicochemical properties but different structure.

> Is synthetic accessibility really important, why?

Yes, as it is a measure of how easy it is to make a compound.

> Is chirality easy to make with _synthetic chemistry_ (not biocatalysis)?

No, as it requires a chiral centre.
Chiral heterocatalysts are not common.
Chirex columns are used to separate enantiomers.
Crystallisation can be used to separate enantiomers —remember scene in Breaking Bad where they smash up crystals...

> Different isomers of racemic compounds bind in different poses within the same crystal. 
> Do these ligands both bind together in the same protein macromolecule at the same time, 
> sequentially in time or on different macromolecules in the crystal lattice?

If the bind the same site, then they are not bound at the same time in a given macromolecule,
they are bound in different macromolecules in the crystal lattice.
In solution they will bind sequentially in time.


### Chose a compound
```python
mol_name: str =  'caffeine'
smiles: str =  'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
```

### 2nd Qeusiton cell

> what is 'RDKit'?

It is a cheminformatics module in Python.

> What is a SMILES (Simplified Molecular-Input Line-entry System) pattern?

A string representation of a molecule.

> What does the function `AllChem.ComputeGasteigerCharges` do?

It computes partial charges on atoms. Marsili-Gasteiger partial charges.
Formal charges are the charges on atoms in a Lewis structure (+1, -1, 0, etc.).

> In addition to Marsili-Gasteiger partial charges, there is another form of partial charges, what is it?


Mulliken Charges: Derived from quantum mechanical calculations, Mulliken charges are obtained by analyzing the electron density in molecular orbitals. They are one of the earliest and simplest methods but can be sensitive to the basis set used in the calculations.

Hirshfeld Charges: These charges are based on a stockholder partitioning of the electron density, offering a more physically realistic distribution of electron density across atoms in a molecule compared to Mulliken charges.

Natural Bond Orbital (NBO) Charges: NBO analysis provides charges based on the concept of "natural" atomic orbitals, which are more localized and can provide a more chemically intuitive picture of charge distribution.

Electrostatic Potential (ESP) Charges: These are derived from fitting the molecular electrostatic potential, which is calculated from quantum mechanical methods. ESP charges are often used in force field development and are considered to be more accurate in representing the electrostatic interactions.

CHELPG (Charges from Electrostatic Potentials using a Grid-based method): Similar to ESP charges, CHELPG charges are derived from the electrostatic potential but use a specific grid-based method for calculation.

RESP (Restrained ElectroStatic Potential) Charges: Developed for use in the AMBER force field, RESP charges are derived from fitting to the electrostatic potential while applying a restraint to keep the charges close to ab initio values.

CM1A/CM2 Charges (Charge Model 1A/2): These are semi-empirical charge models developed by Cramer and Truhlar, designed to produce partial charges that are consistent with dipole moments and polarizabilities.

QTAIM (Quantum Theory of Atoms in Molecules) Charges: These charges are derived from the QTAIM theory, which provides a way of partitioning the electron density in a molecule into atomic basins.

### Third question cell

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

### Fourth question cell

> What is the SMILES of the compound you searched in SmallWorld? And what was the distance to a purchasable analogue?

Chose what form of inappropriateness you want to use:

Fentanyl: CCN(CC)C(=O)C1(CC2=C(C=C1)N(C)C(=O)C2)C3=CC=CC=C3
MDMA: CC(CC1=CC2=C(C=C1)OCO2)NC
Sildenafil: CC1=NN(C)C2=C1C(=O)N(CC3CC4=CC=CC=C4CC3)C(N2C)=O

> Did you find a molecule that has a Wikipedia page but is not in EnamineREAL DB?

Yes, it was a natural products.

> How many small molecules are there?

10^60 they say
Cf. GDB-17: 166 billion molecules up to 17 atoms (https://pubs.acs.org/doi/10.1021/acs.jcim.9b00237)


### Rotation

```python
cen: npt.NDArray[np.float64] = get_centroid(mol.GetConformer())
shift_to_origin: npt.NDArray[np.float64] = create_translation_matrix(*(-cen))
shift_from_origin: npt.NDArray[np.float64] = create_translation_matrix(*(+cen))
rotation: npt.NDArray[np.float64] = create_rotation_matrix([1,0,0], angle=np.pi / 2)
# combine the affine transform matrices by matrix multiplication
_r: npt.NDArray[np.float64] = np.matmul(shift_to_origin, rotation)
rotation_on_spot: npt.NDArray[np.float64] = np.matmul(_r, shift_from_origin)

transform(mol.GetConformer(), rotation_on_spot )
```

### Fifth question cell

> What does a translation look like in a 4x4 matrix?

The last vertical column is the translation vector.

> A dot product between a 3x1 vector and a 4x4 matrix is not possible, what is the missing detail?

A dot product (or matrix multiplication) between a 3x1 vector and a 4x4 matrix is not directly possible under standard matrix multiplication rules because the number of columns in the first matrix (or vector in this case) must match the number of rows in the second matrix. 

This is where the concept of homogeneous coordinates comes into play, particularly in computer graphics and computational geometry.

Homogeneous coordinates are a system of coordinates used in projective geometry, 
as they allow for the representation of points at infinity. They are also extensively used in computer graphics to manage transformations such as translation, rotation, scaling, and perspective projection in a unified way.

Extending Standard Coordinates: A point in 2D Cartesian coordinates (x,y) can be represented in homogeneous coordinates as 
(wx,wy,w), where w is a non-zero scalar known as the weight. For 3D points, a point
(x,y,z) in Cartesian coordinates can be represented as 
(wx,wy,wz,w) in homogeneous coordinates.

Homogenization: In the context of the question, to multiply a 3D vector with a 4x4 matrix (common in 3D transformations), the vector is first converted into a 4D vector in homogeneous coordinates. This is usually done by adding a fourth component, often set to 1, resulting in 
(x,y,z,1). This additional component allows for the inclusion of translation transformations in matrix operations.

Transformation Matrices: In 3D graphics, 4x4 matrices are used to perform transformations. The last row of these matrices, in homogeneous coordinates, often takes the form of
(0,0,0,1), which allows for both linear transformations (like rotation and scaling) and translation.

Converting Back to Cartesian Coordinates: After applying the transformation using the 4x4 matrix, the resulting homogeneous coordinates are often converted back to Cartesian coordinates by dividing each component by the weight (the last component of the vector).

In summary, homogeneous coordinates extend the dimensionality of points by one extra dimension to facilitate easy handling of transformations, particularly those involving translations. This is why they are key in fields like computer graphics and robotics, where such transformations are frequently required.

> How do you rotate a molecule around a point?

Translate to origin, rotate, translate back.

### MCMC Docking

## Understanding MCMC

- **Monte Carlo Markov Chain (MCMC)** is a statistical method for sampling from probability distributions. In ligand-protein docking, it explores possible conformations and orientations of the ligand relative to the protein.
- MCMC generates a chain of samples where each sample depends only on the previous one (Markov property), and the transition probabilities are determined by a Monte Carlo method.

## Initial Setup

- Start with known structures of the target protein and the ligand.
- Define the binding site on the protein, or consider the entire protein surface if the binding site is unknown.

## Random Initialization

- Place the ligand in a random position and orientation within or near the binding site.

## Monte Carlo Sampling

- Perform random moves (translations, rotations, conformational changes) to the ligand.
- After each move, evaluate the energy of the ligand-protein complex using a scoring function that accounts for various molecular interactions.

## Acceptance Criterion

- Use the **Metropolis Criterion** to decide whether to accept or reject each new pose, based on energy changes.
- Lower-energy configurations are more likely to be accepted, but some higher-energy configurations are accepted probabilistically to avoid local minima.

## Convergence and Analysis

- Repeat the sampling process until convergence is reached.
- Analyze the results by examining frequently occurring or low-energy poses, indicative of favorable ligand-protein interactions.

## Refinement and Validation

- Refine the best-scoring poses with more accurate methods.
- Validate the predicted poses using experimental data or other computational methods if available.

MCMC-based docking efficiently explores conformational space and can escape local minima. However, its accuracy depends on the scoring function and the specific implementation of the Monte Carlo algorithm.


# Notebook 2

Both PDBFile (the IO for PDB files) and Modeller (the builder) have a `.topology` and `.positions` attributes.

> What do they look like and what are they describing? (remember `dir` and `type`)

* `.topology` is a graph of the molecule, with atoms and bonds and chains and residues.
* `.positions` is a numpy array of the cartesian coordinates of the atoms inside a OpenMM Quantity like a Pint Quantity.
The latter has a value and a unit. `_value` holds the actual numpy array
Pint module: https://pint.readthedocs.io/en/stable/

### Picking a forcefield

A page on the web says this `forcefield = mma.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')`

But we will use this as we will use implicit water: `forcefield = mma.ForceField('amber14-all.xml', 'implicit/gbn2.xml')`

> What is the difference? And why can't use vacuum? I read on Reddit water is an intersubjective construct...

### Differences in Ligand Binding Energy: Vacuum vs. Solvent

The ligand binding energy in vacuum can differ from that in a solvent due to various physical and chemical interactions. These differences are important in understanding biomolecular interactions.

#### Solvent Effects

- **Polarization**: Polar solvents stabilize interactions through hydrogen bonds and other interactions, absent in vacuum.
- **Dielectric Constant**: The dielectric constant affects electrostatic interactions, which are weaker in high dielectric solvents like water compared to vacuum.

#### Hydrophobic and Hydrophilic Interactions

- **Hydrophobic Effect**: In aqueous solutions, hydrophobic parts of the ligand and protein drive binding by avoiding water, an effect absent in vacuum.
- **Hydrophilic Interactions**: Favorable interactions with the solvent can influence binding affinity, not present in vacuum.

#### Entropy Considerations

- **Solvent Reorganization**: Changes in solvent organization around the ligand and protein contribute to entropy change upon binding.
- **Desolvation**: The displacement of solvent molecules from the binding site during ligand binding contributes to the energetics of binding, absent in vacuum.

#### Specific Solvent Interactions

- Solvents can form specific interactions with the ligand or protein, stabilizing certain conformations not present in vacuum.

#### Computational Modeling Considerations

- **Implicit vs. Explicit Solvent Models**: Different modeling approaches can lead to different predictions about binding energies.
- **Limitations of Vacuum Calculations**: Calculations in vacuum often fail to account for the complex interplay of forces in a physiological environment.

In summary, the interaction with the solvent, screening of electrostatic interactions, conformational changes, 
and solvation/desolvation energetics are key factors in the difference in binding energy between solvent and vacuum environments.

### Error

> What did we do wrong? Take a guess!

We did not add hydrogens to the protein

### Q


> What is this Integrator thing? (cf. http://docs.openmm.org/latest/userguide/theory/04_integrators.html)

An OpenMM Integrator is an algorithm (a symplectic integrator) for integrating a reformulation of Newton's equations of motion.
Newtonian mechanics -> Lagrangian Mechanics -> Langevin Dynamics.
Hamiltonian mechanics (via Velvet and LeapFrog) can also be used.

> What are these forces?


(cf. http://docs.openmm.org/latest/userguide/theory/02_standard_forces.html and theory page in repo)

> 12,000 kcal/mol is a lot of energy. What is going on?

It was not energy minimised

# Lenard-Jones

> We saw the Lenard-Jones potential. What does it do and look like?

https://en.wikipedia.org/wiki/Lennard-Jones_potential
Hockey stick. two terms.

Morse has two terms also, but different. Same shape.
https://en.wikipedia.org/wiki/Morse_potential

> Why is the kinetic energy zero? Even if we _set_ an integrator in Langevin dynamics?

It is at Zero kelvin, local minimum

## Q

> What atom types can you see (cf. https://ambermd.org/antechamber/gaff.html#atomtype or paper for GAFF2)?

Hybridisation state, element symbol, aliphatic or aromatic, apolar or polar, number of hydrogens.
For aromatic nitrogen in an imidazole (cf. histidine) one is protonated and the other is not.

> Why use atom types instead of element symbols?

Geometry is different for different hybridisation states.

> Does an `aromatic sp2 C` form trans or cis isomers?

For a 6 membered ring, the bond angles are 120 degrees, so the hydrogens are cis.
Nobody cares about weird annulenes.

> What about `aliphatic sp2 C`? 

Both

> Is the hydrogen - heavy atom bond the same length?

Shorter by 0.1 Å

# Other questions

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

