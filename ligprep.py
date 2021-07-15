import os, sys

import schrodinger.structure as structure
import schrodinger.pipeline.stages.ligprep as ligprep
import schrodinger.pipeline.pipeio as pipeio


def prepare_ligands(input : str) -> str:
    """
    A calling function that runs the Schrodinger Maestro ligprep on a raw, unprepared ligand file

    Args:
        input (str) : The name of the ligand input file, located in ligands/{input}.sdf

    Returns: 
        output_path (str) : The file path of the prepared ligand after ligprep has finished (located in ligprep dir)
    """

    ligands = os.path.join('ligands', f'{input}.sdf')

    if not os.path.isfile(ligands):
        raise FileNotFoundError(f"The ligand input file {ligands} does not exist!")

    ligands = os.path.abspath(ligands)

    start_dir = os.getcwd()
    work_dir = os.path.join('ligprep', input)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

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
    infile_list = [ligands]
    ligandsobj = pipeio.Structures(infile_list)
    lps.setInput(1, 'INPUT1', ligandsobj)

    # Where to save output
    output_file = f'{input}_prepared.maegz'

    lps.setOutputName(1, f'{input}_prepared')

    # Outputs (dictionary of output objects)
    lps.run()

    prepfile = f'{input}_prepared-001.maegz'
    os.system(f'mv {prepfile} {output_file}')

    os.chdir(start_dir)

    output_path = os.path.join(work_dir, output_file)
    return output_path

if __name__ == '__main__':
    prepare_ligands(sys.argv[1])