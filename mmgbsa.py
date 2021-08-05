import sys, os, argparse, subprocess, shutil

if __name__ == '__main__':
    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger MMGBSA on a directory with a Schrodinger Pose File')

    parser.add_argument('ligands', help='[string] The name of the the set of ligands to run MMGBSA on')

    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)
    parser.add_argument('--do_pocket', help='[bool] Run the MMGBSA on the pocket pose file, rather than the full geometry', default=False, action='store_true')

    args = parser.parse_args(sys.argv[1:])

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    ligands = args.ligands
    ncore = args.ncore
    do_pocket = args.do_pocket

    ## => Run MMGBSA <= ##
    ligands = f'{ligands}_val'

    if do_pocket:
        posefile = os.path.abspath(os.path.join('docking', ligands, 'pocket_pv.maegz'))
        output = 'pocket-out.csv'
    else:
        posefile = os.path.abspath(os.path.join('docking', ligands, 'dockjob_pv.maegz'))
        output = 'dockjob-out.csv'

    start_dir = os.getcwd()
    work_dir = os.path.join('mmgbsa', ligands)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    mmgbsa_job = subprocess.Popen([f'{schrodinger_path}/prime_mmgbsa', posefile, '-out_type', 'PV', '-job_type', 'ENERGY', '-HOST', f'localhost:{ncore}', '-WAIT'])
    mmgbsa_job.wait()

    shutil.move(output, 'output.csv')
    os.chdir(start_dir)

    if mmgbsa_job.returncode != 0:
        raise Exception("The MMGBSA job has failed.")