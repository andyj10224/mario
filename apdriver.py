import os, sys, argparse, subprocess, shutil

dg_path = os.environ.get('APNETDG')
if dg_path is None: raise Exception("Environment variable $APNETDG is not set.")

schrodinger_path = os.environ.get('SCHRODINGER')
if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='Runs an AP-Net-dG prediction on the pocket of a docked protein-ligand system')
    parser.add_argument('ligands', help='[string] The name of the set of ligands to train and validate AP-Net-dG on')

    args = parser.parse_args(sys.argv[1:])
    ligands = args.ligands

    trainname = f'{ligands}_train'
    valname = f'{ligands}_val'
    modelname = ligands

    traindata = os.path.join('docking', trainname, 'pocket_pv.maegz')
    valdata = os.path.join('docking', valname, 'pocket_pv.maegz')

    ## => Write train and validation data to AP-Net directory <== ##
    train_preprocess = subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', traindata, trainname])
    val_preprocess = subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', valdata, valname])
    train_preprocess.wait()
    val_preprocess.wait()

    ## => Train the AP-Net-dG model, and then evaluate the model on the validation data <= ##
    start_dir = os.getcwd()
    os.chdir(dg_path)
    train_job = subprocess.Popen(['python', 'train_model.py', trainname, valname, 'label', '--name', modelname])
    train_job.wait()
    val_job = subprocess.Popen(['python', 'eval_model.py', modelname, valname, 'label'])
    val_job.wait()
    os.chdir(start_dir)

    ## => Save the AP-Net-dG predictions <= ##
    pred_dir = f'apnetdg/{valname}'
    if not os.path.isdir(pred_dir): os.makedirs(pred_dir)
    shutil.copy(f'{dg_path}/datasets/{valname}_{modelname}_preds.csv', f'apnetdg/{valname}/preds.csv')
    shutil.copy(f'{dg_path}/datasets/{valname}_{modelname}_preds.npy', f'apnetdg/{valname}/preds.npy')

    if train_job.returncode != 0 or train_preprocess.returncode != 0 or val_job.returncode != 0 or val_preprocess.returncode != 0:
        raise Exception("The AP-Net-dG job has failed")