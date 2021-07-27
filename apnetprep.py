import os, argparse, sys

import numpy as np
import pandas as pd
import schrodinger.structure as structure
import schrodinger.protein.captermini as captermini

def get_pockets(protein_st, ligand_st, rC=4.0):
    """
    Returns a single, unified structure containing only the protein pocket atoms (cleaned up) 
    and ligand molecule.

    Parameters:
        protein_st (schrodinger.structure.Structure) = the protein structure object
        ligand_st (schrodinger.structure.Structure) = the ligand structure object
        rC (float) = Cutoff distance for a protein pocket

    Returns:
        pocket_st (schrodinger.structure.Structure) = the pocket-ligand structure
    """
    close_residues = []
    pocket_atoms = []
    for res in protein_st.residue:
        for ratom in res.atom:
            found = False
            rA = np.array([ratom.x, ratom.y, ratom.z])
            for latom in ligand_st.atom:
                rB = np.array([latom.x, latom.y, latom.z])
                if np.linalg.norm(rB - rA) <= rC:
                    close_residues.append(res)
                    pocket_atoms.extend(res.getAtomIndices())
                    found = True
                    break
            if found: break

    pocket_st = protein_st.extract(pocket_atoms, copy_props=True)
    pocket_st = captermini.CapTermini(pocket_st, frag_min_atoms=0).outputStructure()

    return pocket_st

def write_data(posefile, dirname):
    """
    Runs AP-Net-dG to get delta G predictions from machine learning 
    (Adds the dG prediction as a column in the .maegz file)

    Parameters:
        posefile (str) : The path to the posefile (located in docking directory)
        dirname (str) : The name of the directory that will contain the AP-Net-dG predictions
    """

    dg_path = os.environ.get('APNETDG')
    if dg_path is None:
        raise Exception("Environment variable $APNETDG is not set.")

    protein_st = None
    ligand_st = []
    for n, st in enumerate(structure.StructureReader(posefile)):
        if (n == 0):
            protein_st = st
        else:
            ligand_st.append(st)

    pockets = []
    for ligand in ligand_st:
        pockets.append(get_pockets(protein_st, ligand))

    data = {
        'RA' : [],
        'RB' : [],
        'ZA' : [],
        'ZB' : [],
        'system' : []
    }

    for p, l in zip(pockets, ligand_st):
        RA = np.zeros((len(p.atom), 3))
        RB = np.zeros((len(l.atom), 3))
        ZA = np.zeros(len(p.atom))
        ZB = np.zeros(len(l.atom))
        system = l.title
        for n, patom in enumerate(p.atom):
            RA[n] = np.array([patom.x, patom.y, patom.z])
            ZA[n] = patom.atomic_number

        for n, latom in enumerate(l.atom):
            RB[n] = np.array([latom.x, latom.y, latom.z])
            ZB[n] = latom.atomic_number
        
        data['RA'].append(RA)
        data['RB'].append(RB)
        data['ZA'].append(ZA)
        data['ZB'].append(ZB)
        data['system'].append(system)
    
    newdir = f'{dg_path}/datasets/{dirname}'
    if not os.path.isdir(newdir): os.makedirs(newdir)
    df = pd.DataFrame(data=data, dtype='object')
    pickle_path = os.path.join(newdir, 'dimers.pkl')
    df.to_pickle(pickle_path)

    """
    name = st.title
    xyz = st.getXYZ()
    natom = st.atom_total
    net_charge = st.formal_charge
    properties = st.property
    gscore = 0.0
    if (n >= 1):
        gscore = properties.get('r_i_glide_gscore')
    """

if __name__ == '__main__':
    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='Preprocesses a pose file to get ready for AP-Net-DG')

    parser.add_argument('posefile', help='[string] The path to the posefile containing the docked protein-ligand geometries')
    parser.add_argument('output_dir', help='[string] The name of the directory to store the AP-Net-dG results, stored as ./apnetdg/{output_dir}')

    args = parser.parse_args(sys.argv[1:])

    posefile = args.posefile
    dirname = args.output_dir

    write_data(posefile, dirname)