__all__ = ['add_mol_in_pose', 'add_mod_cl']

import pyrosetta
from typing import Callable
import types
import os, gzip, re, pathlib, datetime
import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D
from rdkit_to_params import Params
prc: types.ModuleType = pyrosetta.rosetta.core
prn: types.ModuleType = pyrosetta.rosetta.numeric


def add_mol_in_pose(pose: pyrosetta.Pose, mol: Chem.Mol, name: str='LIG', copy_coordinates: bool=True) -> pyrosetta.Pose:
    """
    Add a Chem.Mol to the pose, by generating a "topology" for it, which in Rosetta are called "ResidueType"
    """
    docked: pyrosetta.Pose = pose.clone()
    topo = Params.from_mol(mol, name=name)
    rts: prc.chemical.ResidueTypeSet = topo.add_residuetype(docked)
    lig: prc.conformation.Residue = prc.conformation.ResidueFactory.create_residue( rts.name_map( 'LIG' ) )
    # fix the coordinates...
    if not copy_coordinates:
        get_pos: Callable[[int, ], Point3D] = lig.GetConformer().GetAtomPosition
        for i in range(mol.GetNumAtoms()):
            p: Point3D = get_pos(i)
            atom: Chem.Atom = topo.mol.GetAtomWithIdx(i)
            atomname: str = atom.GetPDBResidueInfo().GetName().strip()
            xyz = prn.xyzVector_double_t(p.x, p.y, p.z)
            lig.set_xyz(lig.atom_index(atomname), xyz)
    # add the ligand to the pose
    docked.append_residue_by_jump(new_rsd=lig, jump_anchor_residue=docked.num_jump()+1)
    return docked


def add_mod_cl(pose: pyrosetta.Pose, gasteiger:float, xyz:prn.xyzVector_double_t):
    cl = Params.loads(f'''NAME CL
    IO_STRING  CL Z
    TYPE LIGAND
    PROPERTIES CHARGED
    AA UNK
    ATOM CL   Cl  X   {gasteiger:.2f}
    ATOM  V1  VIRT  VIRT    0.00
    ATOM  V2  VIRT  VIRT    0.00
    BOND CL    V1
    BOND CL    V2
    BOND  V1   V2
    CHARGE CL FORMAL -1
    NBR_ATOM CL
    NBR_RADIUS 0.01
    ICOOR_INTERNAL   CL      0.000000    0.000000    0.000000  CL     V1    V2
    ICOOR_INTERNAL    V1     0.000000  180.000000    1.493765  CL     V1    V2
    ICOOR_INTERNAL    V2     0.000000   62.060678    1.415324   V1   CL     V2 ''')
    rts: prc.chemical.ResidueTypeSet = cl.add_residuetype(pose)
    cl_res: prc.conformation.Residue = prc.conformation.ResidueFactory.create_residue(rts.name_map('CL'))
    # fix the coordinates...
    cl_res.set_xyz(1, xyz)
    # add the ligand to the pose
    pose.append_residue_by_jump(new_rsd=cl_res, jump_anchor_residue=pose.num_jump() + 1)



DumpTrajectoryEnergy = prc.energy_methods.DumpTrajectoryEnergy
#formerly: DumpTrajectoryEnergy = prc.scoring.util_methods.DumpTrajectoryEnergy

def enable_trajectory(scorefxn: pyrosetta.ScoreFunction,
                      prefix:str='pose',
                      stride:int=100,
                      gz:bool=False) -> DumpTrajectoryEnergy:
    """
    ... code-block:: python
        experiment_name = 'test123'
        scorefxn = pyrosetta.get_fa_scorefxn()
        new_traj_energy = dtc.enable_trajectory(scorefxn=scorefxn,
                                                  prefix=experiment_name,
                                                  stride=10_000)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
        relax.apply(pose)
        trajs = get_traj_interations(experiment_name)
        trajs.head()
    """

    dt_st = prc.scoring.ScoreType.dump_trajectory
    scorefxn.set_weight(dt_st, 1)
    # ----- make the traj dump less intense -----
    # settable by cmd line options... maybe better
    old_traj_energy = [k for k in scorefxn.all_methods() if isinstance(k, DumpTrajectoryEnergy)][0]
    scorefxn.all_methods().remove(old_traj_energy)
    emo = pyrosetta.rosetta.core.scoring.methods.EnergyMethodOptions()
    if gz:
        # emo.dump_trajectory_gz(True)
        raise RuntimeError('this just outputs a text file with gz suffix... (?!)')
    emo.dump_trajectory_prefix(prefix)
    emo.dump_trajectory_stride(stride)
    new_traj_energy = DumpTrajectoryEnergy(emo)
    scorefxn.all_methods().append(new_traj_energy)
    return new_traj_energy

def get_traj_interations(prefix:str) -> pd.DataFrame:
    raw_trajs = []
    for file in os.listdir():
        if re.match(prefix, file):
            with  open(file,'r') as fh:
                file_content=fh.read()
            modstamp = pathlib.Path(file).stat().st_mtime
            modtime = datetime.datetime.fromtimestamp(modstamp)
            models = re.findall('MODEL', file_content)
            raw_trajs.append(dict(filename=file,
                               time=modtime,
                               n_models=len(models))
                            )
    return pd.DataFrame(raw_trajs).sort_values('time')
