import sys
import os

import schrodinger.structure as structure
import schrodinger.pipeline.stages.prime as prime
import schrodinger.pipeline.pipeio as pipeio

def run_mmgbsa(posefile, outdir):
    """
    Run MMGBSA on a pose file performed from Glide ligand docking

    Args:
        posefile (str) : The path to the pose file to run MMGBSA on
        outdir (str) : The name of the directory where the output is stored (as mmgbsa/{outdir})

    Returns:
        None
    """

    posefile = os.path.abspath(posefile)
    posebase = os.path.splitext(os.path.split(posefile)[-1])[0]

    start_dir = os.getcwd()
    work_dir = os.path.join('mmgbsa', outdir)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    mms = prime.MMGBSAStage(f"{posebase}_mmgbsa")

    mms['USE_MAE_CHARGES'] = False # = boolean(default=False) # Whether to use atom partial charges in input Maestro file in simulations
    mms['OUTPUT_LIG_STRAING'] = False # = boolean(default=False) # Whether to output an estimate of the ligand strain energy
    mms['USE_MEMBRANE'] = False # = boolean(default=False) # Use Prime implicit membrane model in protein and complex simulations

    # Input Ligand Objects
    structobj = pipeio.Structures([posefile])
    mms.setInput(1, 'INPUT1', structobj)

    # Outputs (dictionary of output objects)
    mms.setOutputName(1, f'{posebase}_mmgbsa')
    mms.run()

    os.chdir(start_dir)

if __name__ == '__main__':
    run_mmgbsa(sys.argv[1], sys.argv[2])