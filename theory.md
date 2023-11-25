# Theory

## Ligand–Protein Binding

### Key Concepts

- **Dissociation Constant (KD)**: Represents the affinity of the ligand for the protein, defined as the ratio of the dissociation rate constant (\(k_{off}\)) to the association rate constant (\(k_{on}\)). Lower KD values indicate higher affinity.
- **Logarithmic Relationship with Gibbs Free Energy**: The binding affinity, reflected by KD, has a logarithmic relationship with the Gibbs free energy of binding.
- **IC50**: The concentration of an inhibitor where the binding or biological activity is reduced by half. Used to measure the efficacy of a substance in inhibiting a specific biological or biochemical function.
- **Ligand Binding Mechanisms**:
  - **Lock and Key**: The ligand fits precisely into the binding site of the enzyme or receptor, like a key fits into a lock.
  - **Induced Fit**: The binding site changes shape to accommodate the ligand, suggesting a more dynamic interaction.

## Thermodynamics of Ligand Binding

Ligands bind to proteins when the process is thermodynamically favorable.

### Components

- **Gibbs Free Energy**: The sum of the enthalpy change (\(ΔH\)) minus the product of the entropy change (\(ΔS\)) and temperature (\(T\)). A negative Gibbs free energy (\(ΔG\)) indicates a spontaneous and favorable process.
  \[ ΔG = ΔH - TΔS \]
- **Enthalpy (\(ΔH\))**: Negative enthalpy indicates that the interactions (like hydrogen bonds, ionic interactions, etc.) during ligand binding are energetically favorable.
- **Entropy (\(ΔS\))**: A measure of disorder or randomness in the system. Positive entropy change is favorable.
  - **Desolvation**: The displacement of water molecules from the binding interface, akin to the chelate effect.
  - **Hydrophobic Effect**: Hydrophobic side chains in proteins tend to repel water, which plays a critical role in protein folding and stability.
  - **Rigification**: The reduction in the degrees of rotational freedom upon binding, leading to a decrease in entropy but can contribute to overall binding affinity.


## Geometry and charge

- **Geometry**: The geometry of a molecule is determined by the spatial arrangement of its atoms, influenced by chemical bonds and hybridization states. [Molecular Geometry](https://en.wikipedia.org/wiki/Molecular_geometry)
- **Charge**: Refers to the electrical charge of a molecule or atom, affected by electron distribution and bonding. [Electric Charge](https://en.wikipedia.org/wiki/Electric_charge)

### Atomic Orbitals and Hybridization

- **Atomic Orbitals**: Regions where electrons are likely to be found. In elements like Carbon (C), Nitrogen (N), and Oxygen (O), the 's' and 'p' orbitals are significant. [Atomic Orbital](https://en.wikipedia.org/wiki/Atomic_orbital)
- **Hybridization**: Mixing of atomic orbitals to form new hybrid orbitals, affecting molecular geometry.
  - **sp2 Hybridization**: Involves one 's' orbital and two 'p' orbitals, forming trigonal planar geometry with bond angles of 120°. Associated with double bonds. [sp2 Hybridization](https://en.wikipedia.org/wiki/Orbital_hybridisation#sp2_hybridisation)
  - **sp3 Hybridization**: Involves one 's' orbital and three 'p' orbitals, forming tetrahedral geometry with bond angles of approximately 109°. Common in molecules with single bonds. [sp3 Hybridization](https://en.wikipedia.org/wiki/Orbital_hybridisation#sp3_hybridisation)

### Valence and Formal Charge

- **Valence**: Indicates the number of electrons available for bonding, inferred from the element's position in the periodic table. [Valence (Chemistry)](https://en.wikipedia.org/wiki/Valence_(chemistry))
- **Formal Charge**: An estimate of the distribution of electric charge within a molecule, assuming equal sharing of bonding electrons. [Formal Charge](https://en.wikipedia.org/wiki/Formal_charge)

### Electronegativity and Partial Charge

- **Electronegativity**: A measure of an atom's ability to attract and bond with electrons, leading to partial charges in polar covalent bonds. [Electronegativity](https://en.wikipedia.org/wiki/Electronegativity)
- **Partial Charge**: Arises due to differences in electronegativity between atoms in a chemical bond, resulting in slightly positive and negative charges. [Partial Charge](https://en.wikipedia.org/wiki/Electric_charge#Partial_charge)

## Non-bonded Molecular Interactions

### Hydrogen Bonds

- **Donor and Acceptor**:
  - The atom bonded to the hydrogen is the donor, and the atom with the lone pair is the acceptor.
- **Bond Formation**: Hydrogens are attracted to lone pairs, forming hydrogen bonds.
- **Types**: Includes standard, resonance-assisted, low-barrier, and single well hydrogen bonds.
- **Bond Lengths**:
  - Standard hydrogen bonds have a hydrogen–donor distance of about 1.0 Ångströms (Å), and hydrogen–acceptor distance greater than 1.8 Å.
- **Ideal Bond Angle**: The ideal angle between donor, hydrogen, and acceptor is 180 degrees.
- **Angle Variability**: The angle between the antecedent atom, donor/acceptor, and hydrogen/lone pair (LP) varies based on the hybridization state of the atom, generally weaker for lone pairs. [Hydrogen Bond Wikipedia](https://en.wikipedia.org/wiki/Hydrogen_bond)

### Salt Bridges

- Salt bridges are ionic interactions between oppositely charged residues, often found in proteins and nucleic acids. [Salt Bridge Wikipedia](https://en.wikipedia.org/wiki/Salt_bridge_(protein_and_supramolecular))

### π (Pi) Electron Interactions

- **Types**:
  - **π-π Stacking**: Interaction between aromatic rings.
  - **π-S Bonds**: Involving sulfur atoms and π systems.
  - **π-Cations**: Between a π system and a positively charged group.
- **Orientation Factors**: Driven by dipole moments, resulting from the diffuse nature of electron orbitals. [Pi Electron Wikipedia](https://en.wikipedia.org/wiki/Pi_bond). cf. MM vs. QM

### Hydrophobic Interactions

- Hydrophobic interactions occur between nonpolar molecules or parts of molecules, driving processes like protein folding in aqueous environments. [Hydrophobic Effect Wikipedia](https://en.wikipedia.org/wiki/Hydrophobic_effect)

### Halogen-bond

### Polarisibility

Certain orbitals will shift depending on iteractions.

## Molecular mechanics models

![mm](https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/MM_PEF.png/1920px-MM_PEF.png)

[Molecular mechanics](https://en.wikipedia.org/wiki/Molecular_mechanics) uses classical mechanics to model molecular systems.

The potential energy of all systems in molecular mechanics is calculated using force fields.

A force field in molecular mechanics is a computational method used to estimate the potential energy of a system of atoms or coarse-grained particles. It calculates the energy landscape, from which the acting forces on every particle are derived as a gradient of the potential energy with respect to the particle coordinates. [\[source\]](https://en.wikipedia.org/wiki/Force_field_(chemistry))

### Functional Form of Potential Energy
- The potential energy in molecular mechanics includes bonded and nonbonded terms. The general form for the total energy in an additive force field is:
  \[ E_{\text{total}} = E_{\text{bonded}} + E_{\text{nonbonded}} \]
  where \( E_{\text{bonded}} \) includes bond, angle, and dihedral energy, and \( E_{\text{nonbonded}} \) includes electrostatic and van der Waals energy. [\[source\]](https://en.wikipedia.org/wiki/Force_field_(chemistry))

### Bond Stretching
- Bond stretching is often modeled using a Hooke's law formula:
  \[ E_{\text{bond}} = \frac{k_{ij}}{2}(l_{ij} - l_{0,ij})^2 \]
  where \( k_{ij} \) is the force constant, \( l_{ij} \) is the bond length, and \( l_{0,ij} \) is the equilibrium bond length. [\[source\]](https://en.wikipedia.org/wiki/Force_field_(chemistry))

### Electrostatic Interactions
- Electrostatic interactions are represented by a Coulomb energy:
  \[ E_{\text{Coulomb}} = \frac{1}{4\pi \varepsilon_0} \frac{q_i q_j}{r_{ij}} \]
  where \( q_i \) and \( q_j \) are atomic charges and \( r_{ij} \) is the distance between two atoms. [\[source\]](https://en.wikipedia.org/wiki/Force_field_(chemistry))

### Lennard-Jones Potential
- The Lennard-Jones potential models soft repulsive and attractive interactions (van der Waals forces) and is expressed as:
  \[ V_{\text{LJ}}(r) = 4\varepsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right] \]
  where \( r \) is the distance between particles, \( \varepsilon \) is the depth of the potential well, and \( \sigma \) is the distance at which the particle-particle potential energy is zero. [\[source\]](https://en.wikipedia.org/wiki/Lennard-Jones_potential)

- The \( \frac{1}{r^{12}} \) term models Pauli repulsion at short distances due to overlapping electron orbitals, and the \( \frac{1}{r^{6}} \) term models attraction at long-range interactions (London dispersion force). [\[source\]](https://en.wikipedia.org/wiki/Lennard-Jones_potential)

This is a fitted equation hence the bizarre 12 and 6, and there are variants (e.g. ETEN, Buckingham, etc.) or 
alternative models (e.g. Morse, Mie, etc.).

### Hybrid terms

Hybrid terms (derogatively called fudge factors) 
typically involve a combination of empirical data and mathematical functions
to model a give property and are statistically fitted to experimental data.

For example to accurately model protein backbone behavior, force fields may incorporate hybrid terms that reflect the energy associated with different backbone torsional angles, as informed by the Ramachandran plot. 
In order to represent the "energy cost" of adopting certain φ and ψ angles based on statistical frequency.

### Further afield

The above do not capture all the forces involved in molecular mechanics.
For example, the [CHARMM](https://en.wikipedia.org/wiki/CHARMM) force field includes terms for improper torsions, Urey-Bradley terms, and explicit hydrogen bonds.
[Drude particle](https://en.wikipedia.org/wiki/Drude_particle) is a particle used to model polarisability.

## Beyond Molecular Mechanics

In MM atoms are particles, but in QM they are waves!
Or technically, the orbitals are described by a wave function.
See [QM](https://en.wikipedia.org/wiki/Quantum_mechanics), [DFT](https://en.wikipedia.org/wiki/Density_functional_theory),
[Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation), [Kohn–Sham equations](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations).

For example, above the orientation of aromatic rings was said to be driven by dipole moments.
This is not considered in MM, but is in QM.

### Further reading

* Quintessential organic chemistry textbook: [Clayden, Organic Chemistry](https://solo.bodleian.ox.ac.uk/permalink/44OXF_INST/35n82s/alma991022208526307026) (chapter 4)
* Quintessential compchem textbook: [Cramer, Essentials of Computational Chemistry](https://solo.bodleian.ox.ac.uk/permalink/44OXF_INST/35n82s/alma991012063199707026) (chapter 2)
* Great online tutorials in compchem and cheminformatics: [teachopencadd](https://projects.volkamerlab.org/teachopencadd/)
* Easy MD article: [Molecular Dynamics Simulation for All](https://www.sciencedirect.com/science/article/pii/S0896627318306846)