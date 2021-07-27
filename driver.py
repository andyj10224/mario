import sys, os, argparse, warnings
from subprocess import Popen

if __name__ == '__main__':

    ## ==> Read in the Arguments ##

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to preprocess protein-ligand systems for further analysis')

    parser.add_argument('pdbid', help='[string] The pdbid for the protein, located at pdbs/{pdb}.pdb')
    parser.add_argument('ligands', help='[string] Name of the ligand set, located at ligands/{ligands}.sdf')

    parser.add_argument('--ncore', help='[bool] Number of CPU cores to allocate to the pipeline', default=1)
    parser.add_argument('--retrieve_pdb', help='[bool] Get the pdb online if the file is not in directory?', default=False, action='store_true')
    parser.add_argument('--docking_precision', help='[string] Which level of precision to run the docking job (HTVS, SP, or XP)', default='SP')
    parser.add_argument('--constraint_type', help='[string] The type of constraint to use for the docking (NONE, CORE, or SHAPE)', default='CORE')
    parser.add_argument('--run_mmgbsa', help='[bool] run mmgbsa?', default=False, action='store_true')
    parser.add_argument('--run_apnet', help='[bool] run apnet?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligands = args.ligands
    ncore = int(args.ncore)
    retrieve_pdb = args.retrieve_pdb
    precision = args.docking_precision
    constraint_type = args.constraint_type
    run_mmgbsa = args.run_mmgbsa
    run_apnet = args.run_apnet

    ## ==> Start Pipeline <== ##

    # Schrodinger Environmental Variable (MUST be set)
    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Allows prepwizard and ligprep to be run in parallel
    if retrieve_pdb:
        prepwizard_job = Popen([f'{schrodinger_path}/run', 'prepwizard.py', pdbid, '--retrieve_pdb'])
    else:
        prepwizard_job = Popen([f'{schrodinger_path}/run', 'prepwizard.py', pdbid])

    ligprep_cores = ncore - 1
    if ncore == 1:
        prepwizard_job.wait()
        ligprep_cores = 1
    
    ligprep_job = Popen([f'{schrodinger_path}/run', 'ligprep.py', ligands, '--ncore', str(ligprep_cores)])

    # Wait for prepwizard to finish
    prepwizard_job.wait()

    # After prepwizard job finishes, we can go ahead and make the grid (using files in the prepwizard directory)
    prepared_proteins = []
    prepwizard_dir = os.path.join('prepwizard', pdbid)
    for entry in os.scandir(prepwizard_dir):
        if entry.path.endswith(".mae"):
            prepared_proteins.append(entry.path)

    # Wait for ligprep to finish
    ligprep_job.wait()

    # Get the ligand file
    prepligfile = ""
    ligdir = os.path.join("ligprep", ligands)
    for entry in os.scandir(ligdir):
        if entry.path.endswith(".maegz"):
            prepligfile = entry.path

    # Run grid jobs (for every protein geometry we get)
    for protein_file in prepared_proteins:
        Popen([f'{schrodinger_path}/run', 'gridgen.py', protein_file, pdbid, '--ncore', str(ncore)]).wait()

    # Get the grid files
    grid_files = []
    grid_dir = os.path.join('grids', pdbid)
    for entry in os.scandir(grid_dir):
        if entry.path.endswith(".zip"):
            grid_files.append(entry.path)

    # Start docking jobs
    posedirs = []
    for gridfile in grid_files:
        basename = os.path.splitext(os.path.split(gridfile)[-1])[0]
        dock_dirname = f'{pdbid}_{basename}_{ligands}'
        posedirs.append(os.path.join('docking', dock_dirname))
        refligand = os.path.join('prepwizard', pdbid, f'{basename}_ligand.mae')
        Popen([f'{schrodinger_path}/run', 'dock.py', gridfile, prepligfile, dock_dirname, '--ncore', str(ncore), \
             '--refligand', refligand, '--constraint_type', constraint_type, '--precision', precision]).wait()

    # Read the output of the docking jobs
    pose_files = []
    for pose_dir in posedirs:
        for entry in os.scandir(pose_dir):
            if entry.path.endswith('.maegz'):
                pose_files.append(entry.path)
    
    if run_mmgbsa:
        for posefile in pose_files:
            basename = os.path.splitext(os.path.split(gridfile)[-1])[0]
            Popen([f'{schrodinger_path}/run', 'mmgbsa.py', posefile, basename, '--ncore', str(ncore)]).wait()

    if run_apnet:
        warnings.warn('AP-Net-dG has not been interfaced with the pipeline yet. Code exiting.')