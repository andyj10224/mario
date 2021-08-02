import sys, os, argparse, subprocess

if __name__ == '__main__':
    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger MMGBSA on a directory with a Schrodinger Pose File')

    parser.add_argument('posedir', help='[string] Name of the directory where the posefile is stored, (./docking/{posedir})')
    parser.add_argument('outdir', help='[string] Name of the directory to store MMGBSA outputs (saved as ./mmgbsa/{outdir}')

    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)
    parser.add_argument('--do_pocket', help='[bool] Run the MMGBSA on the pocket pose file, rather than the full geometry', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    posedir = args.posedir
    outdir = args.outdir
    ncore = args.ncore
    do_pocket = args.do_pocket

    ## => Run MMGBSA <= ##

    if do_pocket:
        posefile = os.path.abspath(os.path.join('docking', posedir, 'pocket_pv.maegz'))
    else:
        posefile = os.path.abspath(os.path.join('docking', posedir, 'dockjob_pv.maegz'))

    start_dir = os.getcwd()
    work_dir = os.path.join('mmgbsa', outdir)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    subprocess.Popen([f'{schrodinger_path}/prime_mmgbsa', posefile, '-out_type', 'PV', '-job_type', 'ENERGY', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()
    os.chdir(start_dir)