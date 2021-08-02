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

    args = parser.parse_args(sys.argv[1:])
    trainname = args.trainname
    valname = args.valname

    traindata = os.path.join('docking', trainname, 'pocket_pv.maegz')
    valdata = os.path.join('docking', valname, 'pocket_pv.maegz')

    ## => Write train and validation data to AP-Net directory <== ##
    subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', traindata, trainname]).wait()
    subprocess.Popen([f'{schrodinger_path}/run', 'helper/aputil.py', valdata, valname]).wait()

    ## => Train the AP-Net-dG predictions on the model <= ##
    start_dir = os.getcwd()
    os.chdir(dg_path)
    subprocess.Popen(['python', 'train_model.py', trainname, valname, 'label']).wait()
    os.chdir(start_dir)

    ## => Write the AP-Net-dG predictions back to the posefiles <= ##