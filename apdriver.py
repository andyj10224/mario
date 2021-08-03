import os, sys, argparse, subprocess, shutil

dg_path = os.environ.get('APNETDG')
if dg_path is None: raise Exception("Environment variable $APNETDG is not set.")

schrodinger_path = os.environ.get('SCHRODINGER')
if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##
    parser = argparse.ArgumentParser(description='Runs an AP-Net-dG prediction on the pocket of a docked protein-ligand system')

    parser.add_argument('trainname', help='[string] The name of the directory containing training data, located in docking dir')
    parser.add_argument('valname', help='[string] The name of the directory containing validation data, located in docking dir')
    parser.add_argument('modelname', help='[string] The name of the directory to save the model, saved at {dg_path}/models/{modelname}/model_best.h5')

    args = parser.parse_args(sys.argv[1:])
    trainname = args.trainname
    valname = args.valname
    modelname = args.modelname

    traindata = os.path.join('docking', trainname, 'pocket_pv.maegz')
    valdata = os.path.join('docking', valname, 'pocket_pv.maegz')

    ## => Write train and validation data to AP-Net directory <== ##
    subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', traindata, trainname]).wait()
    subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', valdata, valname]).wait()

    ## => Train the AP-Net-dG model, and then evaluate the model on the validation data <= ##
    start_dir = os.getcwd()
    os.chdir(dg_path)
    subprocess.Popen(['python', 'train_model.py', trainname, valname, 'label', '--name', modelname]).wait()
    subprocess.Popen(['python', 'eval_model.py', modelname, valname, 'label']).wait()
    os.chdir(start_dir)

    ## => Save the AP-Net-dG predictions <= ##
    pred_dir = f'apnetdg/{valname}'
    if not os.path.isdir(pred_dir): os.makedirs(pred_dir)
    shutil.copy(f'{dg_path}/datasets/{valname}_{modelname}_preds.csv', f'apnetdg/{valname}/preds.csv')
    shutil.copy(f'{dg_path}/datasets/{valname}_{modelname}_preds.npy', f'apnetdg/{valname}/preds.npy')