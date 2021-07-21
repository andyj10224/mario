import sys, os, argparse, warnings
from subprocess import Popen

def run_pipeline(pdbid, ligand, retrieve_pdb, ncore, run_mmgbsa, run_apnet):
    """
    The function used to call the pipeline, the central information flow hub.

    Args:
        pdbid (str) : The pdbid of the protein input file, located in pdbs/{pdbid}.pdb
        ligand (str) : The name of the ligand input file, located in ligands/{ligand}.sdf
        retrieve_pdb (bool) : Whether or not to retrieve the pdb if pdbs/{pdbid}.pdb is not found
        ncore (int) : Number of CPU cores to use
        mmgbsa (bool) : Run mmgbsa on all the docked ligand pose files
        apnet (bool) : Run apnet on all the docked ligand pose files

    Return: 
        None
    """

    # Schrodinger Environmental Variable (MUST be set)
    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Allows prepwizard and ligprep to be run in parallel
    prepwizard_job = Popen([f'{schrodinger_path}/run', 'prepwizard.py', pdbid, str(retrieve_pdb)])

    if ncore == 1: prepwizard_job.wait()
    
    ligprep_job = Popen([f'{schrodinger_path}/run', 'ligprep.py', ligand])

    prepwizard_job.wait()

    # After prepwizard job finishes, we can go ahead and make the grid (using files in the prepwizard directory)
    prepared_proteins = []
    prepwizard_dir = os.path.join('prepwizard', pdbid)
    for entry in os.scandir(prepwizard_dir):
        if entry.path.endswith(".mae"):
            prepared_proteins.append(entry.path)

    if ncore == 1: ligprep_job.wait()

    # Run grid jobs independently (for every protein geometry we get)
    grid_processes = []
    count = 0
    for protein_file in prepared_proteins:
        gp = Popen([f'{schrodinger_path}/run', 'gridmaker.py', protein_file, pdbid])
        grid_processes.append(gp)
        count += 1
        if count >= ncore - 1:
            gp.wait()

    # Wait for ligprep and gridprep to finish before docking
    ligprep_job.wait()

    # Get the ligand file
    prepligfile = ""
    ligdir = os.path.join("ligprep", ligand)
    for entry in os.scandir(ligdir):
        if entry.path.endswith(".maegz"):
            prepligfile = entry.path

    for gp in grid_processes:
        gp.wait()

    # Get the grid files
    grid_files = []
    grid_dir = os.path.join('grids', pdbid)
    for entry in os.scandir(grid_dir):
        if entry.path.endswith(".zip"):
            grid_files.append(entry.path)

    # Start docking jobs
    dock_processes = []
    count = 0
    for gridfile in grid_files:
        dp = Popen([f'{schrodinger_path}/run', 'docker.py', gridfile, prepligfile, f'{pdbid}_{ligand}'])
        dock_processes.append(dp)
        count += 1
        if count >= ncore:
            dp.wait()

    # Make sure the jobs finish
    for dp in dock_processes:
        dp.wait()

    # Make sure docking jobs finish
    pose_files = []
    pose_dir = os.path.join('docking', f'{pdbid}_{ligand}')
    for entry in os.scandir(pose_dir):
        if entry.path.endswith('.maegz'):
            pose_files.append(entry.path)
    
    if run_mmgbsa:
        count = 0
        for posefile in pose_files:
            mmp = Popen(['f{schrodinger_path}/run', 'mmgbsa.py', posefile, f'{pdbid}_{ligand}'])
            count += 1
            if count >= ncore:
                mmp.wait()

    if run_apnet:
        warnings.warn('AP-Net-dG has not been interfaced with the pipeline yet. Code exiting.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MARIO: A pipeline to preprocess protein-ligand systems for further analysis')

    parser.add_argument('pdbid', help='[string] The pdbid for the crystallographic protein-ligand geometry on Protein Data Bank')
    parser.add_argument('ligfile', help='[string] Path to the ligand file containing raw data')

    parser.add_argument('--ncore', help='[bool] Number of CPU cores to allocate to the pipeline', default=1)
    parser.add_argument('--retrieve_pdb', help='[bool] Get the pdb online if the file is not in directory?', default=False, action='store_true')
    parser.add_argument('--run_mmgbsa', help='[bool] run mmgbsa?', default=False, action='store_true')
    parser.add_argument('--run_apnet', help='[bool] run apnet?', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    pdbid = args.pdbid
    ligfile = args.ligfile
    ncore = int(args.ncore)
    retrieve_pdb = args.retrieve_pdb
    run_mmgbsa = args.run_mmgbsa
    run_apnet = args.run_apnet

    run_pipeline(pdbid, ligfile, retrieve_pdb, ncore, run_mmgbsa, run_apnet)