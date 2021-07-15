import sys, os, argparse, warnings

import prepwizard as pw
import ligprep as lp
import gridmaker as gr
import docker as dk
import apnetprep as ap
import mmgbsa

def run_pipeline(pdbid : str, liginput : str, retrieve_pdb : bool, run_mmgbsa : bool, run_apnet : bool):
    """
    The function used to call the pipeline, the central information flow hub.

    Args:
        pdbid (str) : The pdbid of the protein input file, located in pdbs/{pdbid}.pdb
        ligand (str) : The name of the ligand input file, located in ligands/{ligand}.sdf
        retrieve_pdb (bool) : Whether or not to retrieve the pdb if pdbs/{pdbid}.pdb is not found
        mmgbsa (bool) : Run mmgbsa on all the docked ligand pose files
        apnet (bool) : Run apnet on all the docked ligand pose files

    Return: 
        None
    """

    prepared_proteins = pw.prepare_pdb(pdbid, retrieve_pdb)
    prepared_ligands = lp.prepare_ligands(liginput)


    grid_files = []
    for protein_file in prepared_proteins:
        grid_files = gr.make_grid(protein_file, pdbid)
        grid_files.extend(grid_files)

    pose_files = []
    for gridfile in grid_files:
        posefile = dk.dock(gridfile, prepared_ligands, f'{pdbid}_{liginput}')
        pose_files.append(posefile)
    
    if run_mmgbsa:
        for posefile in pose_files:
            mmgbsa.run_mmgbsa(posefile, f'{pdbid}_{liginput}')

    if run_apnet:
        warnings.warn('AP-Net-dG has not been interfaced with the pipeline yet. Code exiting.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to preprocess protein-ligand systems for further analysis')

    parser.add_argument('pdbid', help='[string] The pdbid for the crystallographic protein-ligand geometry on Protein Data Bank')
    parser.add_argument('ligfile', help='[string] Path to the ligand file containing raw data')

    parser.add_argument('--retrieve_pdb', help='[bool] Get the pdb online if the file is not in directory?', default=False, action='store_true')
    parser.add_argument('--run_mmgbsa', help='[bool] run mmgbsa?', default=False, action='store_true')
    parser.add_argument('--run_apnet', help='[bool] run apnet?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligfile = args.ligfile
    retrieve_pdb = args.retrieve_pdb
    run_mmgbsa = args.run_mmgbsa
    run_apnet = args.run_apnet

    run_pipeline(pdbid, ligfile, retrieve_pdb, run_mmgbsa, run_apnet)