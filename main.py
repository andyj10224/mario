import sys

import os

import prepwizard as pw
import ligandremover as lr
import ligprep as lp
import gridmaker as gr
import docker as dk

def step1(pdbid : str, lig_input : str) -> tuple:

    prepared_pdbs = pw.prepare_pdb(pdbid)
    lp.prepare_ligands(lig_input)
    prep_ligfile = "CLEANED_LIGANDS-001.maegz"

    grid_files = []
    for pdbfile in prepared_pdbs:
        ligand_free_file = lr.remove_ligands(pdbfile)
        files = gr.make_grid(pdbfile, ligand_free_file, pdbid)
        grid_files.extend(files)

    return (grid_files, prep_ligfile)

def step2(grid_files : list, prep_ligfile : str) -> None:
    for gridfile in grid_files:
        dk.dock(gridfile, prep_ligfile)

def finalize(pdbid : str):

    if not os.path.isdir(pdbid):
        os.system(f'mkdir {pdbid}')

    if not os.path.isdir(f'{pdbid}/logs'):
        os.system(f'mkdir {pdbid}/logs')
    
    os.system(f'mv *.log {pdbid}/logs')

    if not os.path.isdir(f'{pdbid}/inputs'):
        os.system(f'mkdir {pdbid}/inputs')
    
    os.system(f'mv *.inp {pdbid}/inputs')

    if not os.path.isdir(f'{pdbid}/grids'):
        os.system(f'mkdir {pdbid}/grids')

    os.system(f'mv *.zip {pdbid}/grids')

    if not os.path.isdir(f'{pdbid}/outputs'):
        os.system(f'mkdir {pdbid}/outputs')
    
    os.system(f'mv *.mae {pdbid}/outputs')
    os.system(f'mv *.maegz {pdbid}/outputs')

if __name__ == '__main__':
    grid_files, prep_ligfile = step1(sys.argv[1], sys.argv[2])
    step2(grid_files, prep_ligfile)
    finalize(sys.argv[1])