import sys
import os

import schrodinger.structure as structure
import schrodinger.pipeline.stages.prime as prime
import schrodinger.pipeline.pipeio as pipeio

def run_mmgbsa(pdbid : str, infile : str):
    mms = prime.MMGBSAStage("MMGBSA")

    mms['USE_MAE_CHARGES'] = False # = boolean(default=False) # Whether to use atom partial charges in input Maestro file in simulations
    mms['OUTPUT_LIG_STRAING'] = False # = boolean(default=False) # Whether to output an estimate of the ligand strain energy
    mms['USE_MEMBRANE'] = False # = boolean(default=False) # Use Prime implicit membrane model in protein and complex simulations

    # Input Ligand Objects
    infile_list = [infile]
    structobj = pipeio.Structures(infile_list)
    mms.setInput(1, 'INPUT1', structobj)

    # Outputs (dictionary of output objects)
    mms.setOutputName(1, 'CLEANED_OUTPUTS')
    outputs = mms.run()

    mmgbsa = os.path.join(pdbid, 'mmgbsa')

    if not os.path.isdir(mmgbsa):
        os.system(f'mkdir {mmgbsa}')

    os.system(f'mv *.maegz {mmgbsa}')
    os.system(f'mv *.log {mmgbsa}')
    os.system(f'mv *.mae {mmgbsa}')
    os.system(f'mv *.csv {mmgbsa}')


if __name__ == '__main__':
    run_mmgbsa(sys.argv[1], sys.argv[2])