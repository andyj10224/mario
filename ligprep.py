import os, sys, argparse, subprocess, warnings, math
import schrodinger.structure as structure

datatypes = ['s', 'i', 'r']
labels = ['Ki', 'IC50', 'Kd', 'EC50']

label_keys = [f'{d}_sd_{l}_(nM)' for d in datatypes for l in labels]

def convert_to_float(data):
    """
    Converts a data point to a float

    data (int, float, or str) : A piece of data from the sdf describing the activity of the ligand

    Returns:
        If conversion was successful, the processed float value.
        Otherwise, None
    """

    if isinstance(data, float):
        return data
    elif isinstance(data, int):
        return float(data)
    elif isinstance(data, str):
        numbers = []
        curr_number = ''
        for char in data:
            if char.isnumeric() or char == '.':
                curr_number += char
            else:
                if curr_number != '':
                    numbers.append(float(curr_number))
                    curr_number = ''

        if len(numbers) == 1:
            return numbers[0]
        else:
            return None
    else:
        return None

def train_val_split(ligname):
    """
    Splits a ligand, stored as ligands/{ligname}.sdf, into training and validation sets

    Parameters:
        ligname (str) : The name of the set of ligands, stored as ligands/{ligname}.sdf
    """
    ligandfile = f'ligands/{ligname}.sdf'
    ligands = []
    ligreader = structure.StructureReader(ligandfile)
    breakout = False
    to_continue = False
    while (True):
        try:
            lig = next(ligreader)
        except Exception as exception:
            if isinstance(exception, StopIteration):
                breakout = True
            else:
                to_continue = True

        if breakout: break
        
        if to_continue:
            to_continue = False
            continue

        label = None
        for k in label_keys:
            if k in lig.property.keys():
                if label is None:
                    label = convert_to_float(lig.property[k])
                else:
                    test_label = convert_to_float(lig.property[k])
                    if test_label is not None:
                        label = min(label, test_label)
        if label is None:
            warnings.warn("Not all of the data has a trainable label.")
            continue

        label = 9.0 - math.log(label, 10.0)
        lig.property['r_sd_training_label'] = label
        ligands.append(lig)
    
    ligands = sorted(ligands, key= lambda x: x.property['r_sd_training_label'])

    train_writer = structure.StructureWriter(f'ligands/{ligname}_train.sdf')
    val_writer = structure.StructureWriter(f'ligands/{ligname}_val.sdf')

    for n, lig in enumerate(ligands):
        if n % 5 == 1:
            val_writer.append(lig)
        else:
            train_writer.append(lig)
    
    train_writer.close()
    val_writer.close()

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger LigPrep Job on a Set of Ligands')

    parser.add_argument('ligands', help='[string] The name of the set of ligands, located at ligands/{ligands}.pdb')
    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])
    ligname = args.ligands
    ncore = int(args.ncore)

    ## ==> Check to make sure everything is there <== ##
    ligands = os.path.join('ligands', f'{ligname}.sdf')
    if not os.path.isfile(ligands): raise FileNotFoundError(f"The ligand input file {ligands} does not exist!")

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    ## ==> Split the ligands into training/validation <== ##
    train_val_split(ligname)

    ## ==> Run the Ligprep Job <== ##
    train_ligands = os.path.join('ligands', f'{ligname}_train.sdf')
    val_ligands = os.path.join('ligands', f'{ligname}_val.sdf')

    train_ligands = os.path.abspath(train_ligands)
    val_ligands = os.path.abspath(val_ligands)

    start_dir = os.getcwd()
    work_dir = os.path.join('ligprep', ligname)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    # Where to save output
    train_output = 'train_prepared.maegz'
    val_output = 'val_prepared.maegz'

    # Run ligprep on train and test data
    subprocess.Popen([f'{schrodinger_path}/ligprep', '-isd', train_ligands, '-omae', train_output, '-epik', '-ph', '7.4', '-pht', '0.1', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()
    subprocess.Popen([f'{schrodinger_path}/ligprep', '-isd', val_ligands, '-omae', val_output, '-epik', '-ph', '7.4', '-pht', '0.1', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()
    
    os.chdir(start_dir)