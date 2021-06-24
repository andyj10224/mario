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

    logs = os.path.join(pdbid, 'logs')
    inputs = os.path.join(pdbid, 'inputs')
    grids = os.path.join(pdbid, 'grids')
    outputs = os.path.join(pdbid, 'outputs')

    if not os.path.isdir(logs):
        os.system(f'mkdir {logs}')
    
    os.system(f'mv *.log {logs}')

    if not os.path.isdir(inputs):
        os.system(f'mkdir {inputs}')
    
    os.system(f'mv *.inp {inputs}')

    if not os.path.isdir(grids):
        os.system(f'mkdir {grids}')

    os.system(f'mv *.zip {grids}')

    if not os.path.isdir(outputs):
        os.system(f'mkdir {outputs}')
    
    os.system(f'mv *.mae {outputs}')
    os.system(f'mv *.maegz {outputs}')

if __name__ == '__main__':
    grid_files, prep_ligfile = step1(sys.argv[1], sys.argv[2])
    step2(grid_files, prep_ligfile)
    finalize(sys.argv[1])