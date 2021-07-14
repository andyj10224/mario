import sys, os, argparse

import prepwizard as pw
import ligandremover as lr
import ligprep as lp
import gridmaker as gr
import docker as dk
import apnetprep as ap

def run_pipeline(pdbid : str, liginput : str, retrieve_pdb : bool) -> tuple:

    if not liginput.endswith('.sdf'):
        raise Exception("Ligand input must be in .sdf format")
    if not os.path.isfile(liginput):
        raise Exception(f"File {liginput} does not exist")

    pdir = os.path.join('proteins', pdbid)
    ldir = os.path.join('ligands', liginput[:-4])
    
    if not os.path.isdir(ldir): os.makedirs(ldir)
    if not os.path.isdir(pdir): os.makedirs(pdir)

    prepared_pdbs = pw.prepare_pdb(pdbid, retrieve_pdb)
    lp.prepare_ligands(liginput)

    # The prepared ligand file
    prep_ligfile = os.path.join('ligands', liginput[:-4], f'{liginput[:-4]}_prepared.maegz')

    grid_files = []
    for pdbfile in prepared_pdbs:
        ligand_free_file = lr.remove_ligands(pdbfile)
        files = gr.make_grid(pdbfile, ligand_free_file, pdbid)
        grid_files.extend(files)

    for gridfile in grid_files:
        dk.dock(gridfile, prep_ligfile, pdbid, liginput)

def prep_apnet():

    return True

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to preprocess protein-ligand systems for further analysis')

    parser.add_argument('pdbid', help='[string] The pdbid for the crystallographic protein-ligand geometry on Protein Data Bank')
    parser.add_argument('ligfile', help='[string] Path to the ligand file containing raw data')

    parser.add_argument('--mmgbsa', help='[bool] run mmgbsa?', default=False)
    parser.add_argument('--retrieve_pdb', help='[bool] Get the pdb online if the file is not in directory?', default=False)

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligfile = args.ligfile
    mmgbsa = args.mmgbsa
    retrieve_pdb = args.retrieve_pdb

    run_pipeline(pdbid, ligfile, retrieve_pdb)