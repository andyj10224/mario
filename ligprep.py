import os, sys

import schrodinger.structure as structure
import schrodinger.pipeline.stages.ligprep as ligprep
import schrodinger.pipeline.pipeio as pipeio


def prepare_ligands(infile : str) -> None:
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
    infile_list = [infile]
    ligandsobj = pipeio.Structures(infile_list)
    lps.setInput(1, 'INPUT1', ligandsobj)

    # Where to save output
    lps.setOutputName(1, f'{infile[:-4]}_prepared')

    # Outputs (dictionary of output objects)
    lps.run()

    prepfile = f'{infile[:-4]}_prepared-001.maegz'
    infile_final = os.path.join('ligands', infile[:-4], f'{infile[:-4]}_raw.sdf')
    prepfile_final = os.path.join('ligands', infile[:-4], f'{infile[:-4]}_prepared.maegz')

    os.system(f'cp {infile} {infile_final}')
    os.system(f'mv {prepfile} {prepfile_final}')