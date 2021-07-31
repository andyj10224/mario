import os, sys, argparse, subprocess, shutil

def eval_data(dirname):

    dg_path = os.environ.get('APNETDG')
    if dg_path is None:
        raise Exception("Environment variable $APNETDG is not set.")
    
    start_dir = os.getcwd()
    store_dir = os.path.join(start_dir, 'apnetdg', dirname)
    if not os.path.isdir(store_dir): os.makedirs(store_dir)

    os.chdir(dg_path)
    subprocess.Popen(['python', 'eval_model.py', 'pocket_model_acsf', dirname, 'label']).wait()
    os.chdir(start_dir)

    # Copy csv and npy AP-Net-dG predictions
    shutil.copy(f'{dg_path}/datasets/{dirname}_pocket_model_acsf_preds.csv', f'{start_dir}/apnetdg/{dirname}/pocket_model_acsf_preds.csv')
    shutil.copy(f'{dg_path}/datasets/{dirname}_pocket_model_acsf_preds.npy', f'{start_dir}/apnetdg/{dirname}/pocket_model_acsf_preds.npy')


if __name__ == '__main__':

    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='Runs an AP-Net-dG prediction on the pocket of a docked protein-ligand system')

    parser.add_argument('posefile', help='[string] The path to the posefile containing the docked protein-ligand geometries')
    parser.add_argument('output_dir', help='[string] The name of the directory to store the AP-Net-dG results, stored as ./apnetdg/{output_dir}')

    args = parser.parse_args(sys.argv[1:])

    posefile = args.posefile
    dirname = args.output_dir

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None:
        raise Exception("Environment variable $SCHRODINGER is not set.")

    ## => Preprocess the data <== ##
    subprocess.Popen([f'{schrodinger_path}/run', 'apnetprep.py', posefile, dirname]).wait()

    ## => Run the AP-Net-dG predictions <= ##
    eval_data(dirname)