import sys, os, argparse, time
from subprocess import Popen

def fail():
    """
    A helper function that stops the pipeline if there is a failure in the middle of the pipeline.
    """
    raise Exception("MAMA MIA!!! MARIO encountered an error.")

def start_timer(name):
    """
    A helper function that returns the start time of a job.

    Parameters:
        name (str) : The name of the stage of the pipeline.

    Returns:
        The starting time of the job.
    """
    print(f'\tStarting {name} Job...\n')
    return time.time()

def end_timer(name, start_time):
    """
    A helper function that returns the elapsed time of a job.

    Parameters:
        name (str) : The name of the stage of the pipeline.
        start_time (float) : The starting time of the job.
    """
    elapsed_time = time.time() - start_time
    print(f'\t{name} Job finished...\n\tTime: {elapsed_time:.2f} seconds\n')

if __name__ == '__main__':

    ## ==> Read in the Arguments ##

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to automate the preprocessing and analysis of protein-ligand systems.')

    parser.add_argument('pdbid', help='[string] The pdbid for the protein, located at pdbs/{pdb}.pdb')
    parser.add_argument('ligid', help='[string] The id of the ligand at the binding site of the protein')
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
    parser.add_argument('--skip_mmgbsa', help='[bool] skip the MMGBSA calculation?', default=False, action='store_true')
    parser.add_argument('--skip_apnet', help='[bool] skip AP-Net-dG training and validation?', default=False, action='store_true')
    parser.add_argument('--skip_analysis', help='[bool] skip running the analysis after AP-Net-dG and MMGBSA?', default=False, action='store_true')
    parser.add_argument('--mmgbsa_do_pocket', help='[bool] run mmgbsa on the binding pocket of the poses only?', default=False, action='store_true')
    parser.add_argument('--mmgbsa_skip_analysis', help='[bool] skip mmgbsa when running the analysis?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligid = args.ligid.upper()
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
    do_mmgbsa = not args.skip_mmgbsa
    do_apnet = not args.skip_apnet
    do_analysis = not args.skip_analysis
    do_pocket = args.mmgbsa_do_pocket
    do_mmgbsa_analysis = not args.mmgbsa_skip_analysis

    ## ==> Start Pipeline <== ##

    # Schrodinger Environmental Variable (MUST be set)
    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    print(f'\tHere we go!!!\n')

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
        if prepwizard_job.returncode != 0: fail()
        end_timer('Prepwizard', prepwizard_start)
        ligprep_cores = 1
    
    if do_ligprep:
        ligprep_start = start_timer('Ligprep')
        ligprep_job = Popen([f'{schrodinger_path}/run', 'ligprep.py', ligands, '--ncore', str(ligprep_cores)])

        # Wait for ligprep to finish
        ligprep_job.wait()
        if ligprep_job.returncode != 0: fail()
        end_timer('Ligprep', ligprep_start)

    # Wait for prepwizard to finish
    if ncore != 1 and do_prepwizard:
        prepwizard_job.wait()
        if prepwizard_job.returncode != 0: fail()
        end_timer('Prepwizard', prepwizard_start)

    # Run grid jobs
    if do_gridgen:
        gridgen_start = start_timer('Gridgen')
        gridgen_job = Popen([f'{schrodinger_path}/run', 'gridgen.py', pdbid, ligid, '--ncore', str(ncore)])
        gridgen_job.wait()
        if gridgen_job.returncode != 0: fail()
        end_timer('Gridgen', gridgen_start)

    # Start docking jobs
    if do_docking:
        docking_start = start_timer("Docking")
        train_dock = Popen([f'{schrodinger_path}/run', 'dock.py', pdbid, ligands, '--ncore', str(ncore), '--refligand', ligid, '--constraint_type', \
            constraint_type, '--precision', precision, '--pocket_cutoff', str(pocket_cutoff), '--training'])
        val_dock = Popen([f'{schrodinger_path}/run', 'dock.py', pdbid, ligands, '--ncore', str(ncore), '--refligand', ligid, '--constraint_type', \
            constraint_type, '--precision', precision, '--pocket_cutoff', str(pocket_cutoff)])
        train_dock.wait()
        val_dock.wait()
        if train_dock.returncode != 0 or val_dock.returncode != 0: fail()
        end_timer("Docking", docking_start)
    
    # Run the MMGBSA jobs
    if do_mmgbsa:
        mmgbsa_start = start_timer("MMGBSA")
        if do_pocket:
            mmgbsa_job = Popen([f'{schrodinger_path}/run', 'mmgbsa.py', ligands, '--ncore', str(ncore), '--do_pocket'])
        else:
            mmgbsa_job = Popen([f'{schrodinger_path}/run', 'mmgbsa.py', ligands, '--ncore', str(ncore)])
        mmgbsa_job.wait()
        if mmgbsa_job.returncode != 0: fail()
        end_timer("MMGBSA", mmgbsa_start)

    # Run the AP-Net training job
    if do_apnet:
        apnet_start = start_timer("AP-NET-DG")
        ap_job = Popen(['python', 'apdriver.py', ligands])
        ap_job.wait()
        if ap_job.returncode != 0: fail()
        end_timer("AP-NET-DG", apnet_start)

    # Run the analysis part of the calculation
    if do_analysis:
        analysis_start = start_timer("ANALYSIS")
        if not do_mmgbsa_analysis:
            analysis_job = Popen([f'{schrodinger_path}/run', 'analysis.py', f'{ligands}_val', '--skip_mmgbsa'])
        else:
            analysis_job = Popen([f'{schrodinger_path}/run', 'analysis.py', f'{ligands}_val'])
        analysis_job.wait()
        if analysis_job.returncode != 0: fail()
        end_timer("ANALYSIS", analysis_start)

    print(f'\tPipeline finished!!! Thank you Mario! Your quest is complete! :)')