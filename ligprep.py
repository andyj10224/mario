import os, sys, argparse, subprocess

import schrodinger.structure as structure
import schrodinger.pipeline.stages.ligprep as ligprep
import schrodinger.pipeline.pipeio as pipeio

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger LigPrep Job on a Set of Ligands')

    parser.add_argument('ligands', help='[string] The name of the set of ligands, located at ligands/{ligands}.pdb')
    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])
    input = args.ligands
    ncore = int(args.ncore)

    ## ==> Run the Ligprep Job <== ##

    ligands = os.path.join('ligands', f'{input}.sdf')

    if not os.path.isfile(ligands):
        raise FileNotFoundError(f"The ligand input file {ligands} does not exist!")

    ligands = os.path.abspath(ligands)

    start_dir = os.getcwd()
    work_dir = os.path.join('ligprep', input)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Where to save output
    output_file = f'{input}_prepared.maegz'

    subprocess.Popen([f'{schrodinger_path}/ligprep', '-isd', ligands, '-omae', output_file, '-epik', '-ph', '7.4', '-pht', '0.1', '-HOST', f'localhost:{ncore}', '-WAIT']).wait()
    
    os.chdir(start_dir)

    """
    lps = ligprep.LigPrepStage("Ligprep")

    lps['UNIQUEFIELD'] = "NONE" # Field to identify unique compound by
    lps['RETITLE'] = False # Whether to set OUTCOMPOUNDFIELD
    lps['OUTCOMPOUNDFIELD'] = "s_m_title" # Field where to store the compound codes (if RETITLE is True).
    lps['STEREO_SOURCE'] = "parities" # parities/geometry
    lps['USE_EPIK'] = True # Whether to use Epik instead of Ionizer.
    lps['METAL_BINDING'] = False # Use Epik metal binding mode.
    lps['RETAIN'] = True # Retain input variant for each compound.
    lps['PH'] = 7.4 # Target pH
    lps['PHT'] = 0.1 # pH threshold
    lps['IONIZE'] = True # Whether to ionize (Ionizer; Epik always ionizes).
    lps['GENERATE_TAUTOMERS'] = True # Whether to generate tautomers.
    lps['MAX_TAUTOMERS'] = 8 # Maximum number of tautomers to generate (Ionizer).
    lps['NEUTRALIZE'] = False # Whether to neutralize before expanding states.
    lps['MAX_STATES'] = 16 # Maximum number of states to generate (Epik).
    lps['NUM_STEREOISOMERS'] = 32 # Generate this many stereoisomers per protonation state.
    lps['MAX_STEREOISOMERS'] = None # Keep this many low-energy stereoisomers per protonation state.
    lps['NRINGCONFS'] = 1 # Ring conformers per ligand
    lps['MIXLIGS'] = False
    lps['RECOMBINE'] = True # Whether to recombine input ligand files
    lps['COMBINEOUTS'] = False # Combine output files
    lps['SKIP_BAD_LIGANDS'] = True
    lps['REGULARIZE'] = False # Whether to standardize input structures before preparing
    lps['SKIP_NOUNIQUE_LIGANDS'] = False # Whether to skip ligands that have no unique field.
    lps['TAUT_SPEC_FILE'] = None # Custom tautomerizer spec file to use.
    lps['NORMALIZE'] = False # Whether to normalize input variants
    lps['OUTFORMAT'] = "mae"

    # Input Ligand Objects
    ligandsobj = pipeio.Structures([ligands])
    lps.setInput(1, 'INPUT1', ligandsobj)

    # Write the input file
    input_files = lps.getInput(1).getFiles()
    print(input_files[0])

    # Where to save output
    output_file = f'{input}_prepared.maegz'

    lps.setOutputName(1, f'{input}_prepared')

    # Run the Stage
    lps.run()

    prepfile = f'{input}_prepared-001.maegz'
    os.system(f'mv {prepfile} {output_file}')
    os.chdir(start_dir)
    """