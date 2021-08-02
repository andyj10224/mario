import sys, os, argparse, warnings, time
from subprocess import Popen

def start_timer(name):
    print(f'\tStarting {name} Job...\n')
    return time.time()

def end_timer(name, start_time):
    elapsed_time = time.time() - start_time
    print(f'\t{name} Job finished...\n\tTime: {elapsed_time:.2f} seconds\n')

if __name__ == '__main__':

    ## ==> Read in the Arguments ##

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to automate the preprocessing and analysis of protein-ligand systems.')

    parser.add_argument('pdbid', help='[string] The pdbid for the protein, located at pdbs/{pdb}.pdb')
    parser.add_argument('ligands', help='[string] Name of the ligand set, located at ligands/{ligands}.sdf')

    parser.add_argument('--ncore', help='[bool] Number of CPU cores to allocate to the pipeline', default=1)
    parser.add_argument('--retrieve_pdb', help='[bool] Get the pdb online if the file is not in directory?', default=False, action='store_true')
    parser.add_argument('--docking_precision', help='[string] Which level of precision to run the docking job (HTVS, SP, or XP)', default='SP')
    parser.add_argument('--constraint_type', help='[string] The type of constraint to use for the docking (NONE, CORE, or SHAPE)', default='NONE')
    parser.add_argument('--pocket_cutoff', help='[float] The cutoff distance (Angstroms) from the ligand that defines the binding pocket', default=5.0)
    parser.add_argument('--skip_prepwizard', help='[bool] skip the protein preparation (if it\'s already done)', default=False, action='store_true')
    parser.add_argument('--skip_ligprep', help='[bool] skip the ligand preparation (if it\'s already done)', default=False, action='store_true')
    parser.add_argument('--skip_gridgen', help='[bool] skip the grid generation (if it\'s already done)', default=False, action='store_true')
    parser.add_argument('--skip_docking', help='[bool] skip the ligand docking (if it\'s already done)', default=False, action='store_true')
    parser.add_argument('--run_mmgbsa', help='[bool] run mmgbsa?', default=False, action='store_true')
    parser.add_argument('--run_apnet', help='[bool] run apnet?', default=False, action='store_true')
    parser.add_argument('--run_analysis', help='[bool] run the analysis after AP-Net-dG and MMGBSA?', default=False, action='store_true')
    parser.add_argument('--mmgbsa_do_pocket', help='[bool] run mmgbsa on the binding pocket of the poses only?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligands = args.ligands
    ncore = int(args.ncore)
    retrieve_pdb = args.retrieve_pdb
    precision = args.docking_precision
    constraint_type = args.constraint_type
    pocket_cutoff = args.pocket_cutoff
    do_prepwizard = not args.skip_prepwizard
    do_ligprep = not args.skip_ligprep
    do_gridgen = not args.skip_gridgen
    do_docking = not args.skip_docking
    do_mmgbsa = args.run_mmgbsa
    do_apnet = args.run_apnet
    do_analysis = args.run_analysis
    do_pocket = args.mmgbsa_do_pocket

    ## ==> Start Pipeline <== ##

    # Schrodinger Environmental Variable (MUST be set)
    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Allows prepwizard and ligprep to be run in parallel
    if do_prepwizard:
        prepwizard_start = start_timer('Prepwizard')
        ligprep_cores = ncore - 1
        if retrieve_pdb:
            prepwizard_job = Popen([f'{schrodinger_path}/run', 'prepwizard.py', pdbid, '--retrieve_pdb'])
        else:
            prepwizard_job = Popen([f'{schrodinger_path}/run', 'prepwizard.py', pdbid])
    else:
        ligprep_cores = ncore

    
    if ncore == 1 and do_prepwizard:
        prepwizard_job.wait()
        end_timer('Prepwizard', prepwizard_start)
        ligprep_cores = 1
    
    if do_ligprep:
        ligprep_start = start_timer('Ligprep')
        ligprep_job = Popen([f'{schrodinger_path}/run', 'ligprep.py', ligands, '--ncore', str(ligprep_cores)])

        # Wait for ligprep to finish
        ligprep_job.wait()
        end_timer('Ligprep', ligprep_start)

    # Wait for prepwizard to finish
    if ncore != 1 and do_prepwizard:
        prepwizard_job.wait()
        end_timer('Prepwizard', prepwizard_start)

    # After prepwizard job finishes, we can go ahead and make the grid (using files in the prepwizard directory)
    prepared_proteins = []
    prepwizard_dir = os.path.join('prepwizard', pdbid)
    for entry in os.scandir(prepwizard_dir):
        if entry.path.endswith(".mae"):
            prepared_proteins.append(entry.path)

    # Run grid jobs (for every protein geometry we get)
    if do_gridgen:
        gridgen_start = start_timer('Gridgen')
        for protein_file in prepared_proteins:
            Popen([f'{schrodinger_path}/run', 'gridgen.py', protein_file, pdbid, '--ncore', str(ncore)]).wait()
        end_timer('Gridgen', gridgen_start)

    # Get the grid files
    grid_files = []
    grid_dir = os.path.join('grids', pdbid)
    for entry in os.scandir(grid_dir):
        if entry.path.endswith(".zip"):
            grid_files.append(entry.path)

    # Start docking jobs
    if do_docking:
        docking_start = start_timer("Docking")
        for gridfile in grid_files:
            basename = os.path.splitext(os.path.split(gridfile)[-1])[0]
            refligand = os.path.join('prepwizard', pdbid, f'{basename}_ligand.mae')
            for subset in ['train', 'val']:
                dock_dirname = f'{pdbid}_{basename}_{ligands}_{subset}'
                # posedirs.append(os.path.join('docking', dock_dirname))
                prepligfile = os.path.join('ligprep', ligands, f'{subset}_prepared.maegz')
                Popen([f'{schrodinger_path}/run', 'dock.py', gridfile, prepligfile, dock_dirname, \
                    '--ncore', str(ncore), '--refligand', refligand, '--constraint_type', \
                    constraint_type, '--precision', precision, '--pocket_cutoff', str(pocket_cutoff)]).wait()
        end_timer("Docking", docking_start)

    # Get the necessary information for MMGBSA and AP-Net-dG
    modelnames = []
    for gridfile in grid_files:
        basename = os.path.splitext(os.path.split(gridfile)[-1])[0]
        modelnames.append(f'{pdbid}_{basename}_{ligands}')
    
    if do_mmgbsa:
        mmgbsa_start = start_timer("MMGBSA")
        for modelname in modelnames:
            dirname = f'{modelname}_val'
            if do_pocket:
                Popen([f'{schrodinger_path}/run', 'mmgbsa.py', dirname, dirname, '--ncore', str(ncore), '--do_pocket']).wait()
            else:
                Popen([f'{schrodinger_path}/run', 'mmgbsa.py', dirname, dirname, '--ncore', str(ncore)]).wait()
        end_timer("MMGBSA", mmgbsa_start)

    if do_apnet:
        apnet_start = start_timer("AP-NET-DG")
        for modelname in modelnames:
            Popen(['python', 'apdriver.py', f'{modelname}_train', f'{modelname}_val', modelname]).wait()
        end_timer("AP-NET-DG", apnet_start)

    if do_analysis:
        analysis_start = start_timer("ANALYSIS")
        for modelname in modelnames:
            Popen(['python', 'analysis.py', f'{modelname}_val']).wait()
        end_timer("ANALYSIS", analysis_start)