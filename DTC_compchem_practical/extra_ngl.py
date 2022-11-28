def add_neighbors(view, selection:str, comp_id:int=0, radius:float=5, style:str = 'hyperball',
                  color:str='gainsboro'):
    view._js(f"""const comp = this.stage.compList[{comp_id}]
                 const target_sele = new NGL.Selection('{selection}');
                 const radius = {radius};
                 const neigh_atoms = comp.structure.getAtomSetWithinSelection( target_sele, radius );
                 const resi_atoms = comp.structure.getAtomSetWithinGroup( neigh_atoms );
                 comp.addRepresentation( "{style}", {{sele: resi_atoms.toSeleString(),
                 									  colorValue: "{color}",
                                                      multipleBond: true
                                                                   }});
                 comp.addRepresentation( "contact", {{sele: resi_atoms.toSeleString()}});
            """)
