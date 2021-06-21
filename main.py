import sys

import prepwizard as pw
import ligprep as lp
import gridmaker as gr
import docker as dk

def step1(pdbid : str, lig_input : str) -> tuple:

    prepared_pdbs = pw.prepare_pdb(pdbid)
    lp.prepare_ligands(lig_input)
    prep_ligfile = "OUTPUT-001.maegz"

    grid_files = []
    for pdbfile in prepared_pdbs:
        files = gr.make_grid(pdbfile, pdbid)
        grid_files.extend(files)

    return (grid_files, prep_ligfile)

def step2(grid_files : list, prep_ligfile : str) -> None:
    for gridfile in grid_files:
        dk.dock(gridfile, prep_ligfile)

if __name__ == '__main__':
    grid_files, prep_ligfile = step1(sys.argv[1], sys.argv[2])
    step2(grid_files, prep_ligfile)