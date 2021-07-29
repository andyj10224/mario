import os, sys, argparse, subprocess

import schrodinger.structure as structure
import schrodinger.pipeline.stages.ligprep as ligprep
import schrodinger.pipeline.pipeio as pipeio

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger LigPrep Job on a Set of Ligands')

    parser.add_argument('ligands', help='[string] The name of the set of ligands, located at ligands/{ligands}.pdb')
    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])
    input = args.ligands
    ncore = int(args.ncore)

    ## ==> Run the Ligprep Job <== ##

    ligands = os.path.join('ligands', f'{input}.sdf')

    if not os.path.isfile(ligands):
        raise FileNotFoundError(f"The ligand input file {ligands} does not exist!")

    ligands = os.path.abspath(ligands)

    start_dir = os.getcwd()
    work_dir = os.path.join('ligprep', input)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Where to save output
    output_file = f'{input}_prepared.maegz'

    subprocess.Popen([f'{schrodinger_path}/ligprep', '-isd', ligands, '-omae', output_file, '-epik', '-ph', '7.4', '-pht', '0.1', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()
    
    os.chdir(start_dir)