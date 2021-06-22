import schrodinger.structure as structure
import schrodinger.protein.findhets as findhets

def remove_ligands(filename : str) -> str:

    struct = structure.StructureReader.read(filename)
    ligand_lists = findhets.find_hets(struct, include_metals=False, include_hydrogens=True)

    atoms_to_remove = []

    for ligand in ligand_lists:
        for atom in ligand:
            atoms_to_remove.append(atom)

    struct.deleteAtoms(atoms_to_remove)

    newfile = f'{filename[:-4]}_no_ligand.mae'
    structure.StructureWriter.write(struct, newfile)

    return newfile