import os, argparse, sys

import numpy as np
import pandas as pd
import schrodinger.structure as structure

def write_data(posefile, dirname):
    """
    Writes the data from a posefile to make a AP-Net-dG dataset

    Parameters:
        posefile (str) : The path to the posefile (located in docking directory)
        dirname (str) : The name of the directory to store the AP-Net-dG inputs (stored as $APNETDG/datasets/{dirname}/dimers.pkl)
    """

    ## ==> Make sure everything is there <== ##
    dg_path = os.environ.get('APNETDG')
    if dg_path is None: raise Exception("Environment variable $APNETDG is not set.")

    structures = []
    for n, st in enumerate(structure.StructureReader(posefile)):
        structures.append(st)

    data = {
        'RA' : [],
        'RB' : [],
        'ZA' : [],
        'ZB' : [],
        'label' : [],
        'system' : []
    }

    protein = structures[0]
    ligands = structures[1:]

    for n, ligand in enumerate(ligands):
        RA = np.zeros((len(protein.atom), 3))
        RB = np.zeros((len(ligand.atom), 3))
        ZA = np.zeros(len(protein.atom))
        ZB = np.zeros(len(ligand.atom))
        system = f'{ligand.title}:{n+1}'
        label = ligand.property['r_sd_training_label']

        for n, patom in enumerate(protein.atom):
            RA[n] = np.array([patom.x, patom.y, patom.z])
            ZA[n] = patom.atomic_number

        for n, latom in enumerate(ligand.atom):
            RB[n] = np.array([latom.x, latom.y, latom.z])
            ZB[n] = latom.atomic_number
        
        data['RA'].append(RA)
        data['RB'].append(RB)
        data['ZA'].append(ZA)
        data['ZB'].append(ZB)
        data['label'].append(label)
        data['system'].append(system)

    newdir = f'{dg_path}/datasets/{dirname}'
    if not os.path.isdir(newdir): os.makedirs(newdir)
    df = pd.DataFrame(data=data, dtype='object')
    pickle_path = os.path.join(newdir, 'dimers.pkl')
    df.to_pickle(pickle_path)

if __name__ == '__main__':
    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='Call various utilities to get a posefile ready for AP-Net-dG evaluation')

    parser.add_argument('posefile', help='[string] The path to the posefile containing the docked protein-ligand geometries')
    parser.add_argument('apdir', help='[string] The directory that stores the AP-Net-dG inputs (stored as $APNETDG/datasets/{dirname}/dimers.pkl)')

    args = parser.parse_args(sys.argv[1:])
    posefile = args.posefile
    apdir = args.apdir

    write_data(posefile, apdir)