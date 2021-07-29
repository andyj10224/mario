import sys, os, argparse, subprocess

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

    ## => Run MMGBSA <= ##

    posefile = os.path.abspath(posefile)

    start_dir = os.getcwd()
    work_dir = os.path.join('mmgbsa', outdir)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    subprocess.Popen([f'{schrodinger_path}/prime_mmgbsa', posefile, '-out_type', 'PV', '-job_type', 'ENERGY', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()