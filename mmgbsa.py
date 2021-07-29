import sys, os, argparse, subprocess

import schrodinger.pipeline.stages.prime as prime
import schrodinger.pipeline.pipeio as pipeio

if __name__ == '__main__':
    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger MMGBSA on a Schrodinger Pose File')

    parser.add_argument('posefile', help='[string] The path to the posefile')
    parser.add_argument('output_dir', help='[string] Name of the directory to store MMGBSA outputs (saved as ./mmgbsa/{output_dir}')

    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    posefile = args.posefile
    outdir = args.output_dir
    ncore = args.ncore

    posefile = os.path.abspath(posefile)
    # posebase = os.path.splitext(os.path.split(posefile)[-1])[0]

    start_dir = os.getcwd()
    work_dir = os.path.join('mmgbsa', outdir)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    subprocess.Popen([f'{schrodinger_path}/prime_mmgbsa', posefile, '-out_type', 'PV', '-job_type', 'ENERGY', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()

    """
    mms = prime.MMGBSAStage("MMGBSA")

    mms['USE_MAE_CHARGES'] = False # = boolean(default=False) # Whether to use atom partial charges in input Maestro file in simulations
    mms['OUTPUT_LIG_STRAING'] = False # = boolean(default=False) # Whether to output an estimate of the ligand strain energy
    mms['USE_MEMBRANE'] = False # = boolean(default=False) # Use Prime implicit membrane model in protein and complex simulations

    # Set the MMGBSA number of cores
    mms.setJobOptions(njobs=ncore)

    # Input Ligand Objects
    structobj = pipeio.Structures([posefile])
    mms.setInput(1, 'INPUT1', structobj)

    # Outputs (dictionary of output objects)
    mms.setOutputName(1, "MMGBSA")
    mms.run()

    os.chdir(start_dir)
    """