import subprocess, sys, os, argparse
import helper.pare as pare
import schrodinger.structure as structure
import schrodinger.application.glide.glide as glide

def pocketeer(posefile, pocketfile, rcutoff=4.0):
    """
    Given a posefile, make another posefile containing only the protein pocket (with all ligands)

    Parameters:
        posefile: [str] Name (or path) of the file containing the docked poses
        pocketfile: [str] Name (or path) of the file containing the pocket binding sites
        rcutoff: [float] The cutoff distance from the ligand to define the protein binding pocket
    
    Returns:
        None
    """

    structures = []
    for st in structure.StructureReader(posefile):
        structures.append(st)

    protein = pare.fragment(structures[0], protein_treatment='sidechain')
    ligands = structures[1:]

    merged_ligands = ligands[0].copy()
    for ligand in ligands[1:]:
        merged_ligands.extend(ligand)

    fa_set, h_set = pare.pareProteinWithinDistance(protein, merged_ligands, dmax=rcutoff)
    pocket, renumbering = pare.pareMinimize(structures[0], fa_set, h_set)

    writer = structure.StructureWriter(pocketfile)
    writer.append(pocket)
    writer.extend(ligands)
    writer.close()

if __name__ == '__main__':

     ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger Glide Docking Job on prepared ligands and a grid')

    parser.add_argument('gridfile', help='[string] The path to the gridfile')
    parser.add_argument('ligandfile', help='[string] The path to the ligand set that is being docked')
    parser.add_argument('output_dir', help='[string] The name of the directory to store the Docking job, stored as ./docking/{output_dir}')

    parser.add_argument('--precision', help='[string] Which level of precision to run the docking job (HTVS, SP, or XP)', default='SP')
    parser.add_argument('--refligand', help='[string] The path to the reference ligand (for constraints)', default=None)
    parser.add_argument('--constraint_type', help='[string] The type of constraint to use for the docking (NONE, CORE, or SHAPE)', default='NONE')
    parser.add_argument('--pocket_cutoff', help='[float] The cutoff distance (Angstroms) from the ligand that defines the binding pocket', default=5.0)
    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])

    gridfile = args.gridfile
    ligandfile = args.ligandfile
    outdir = args.output_dir

    precision = args.precision.upper()
    refligand = args.refligand
    constraint_type = args.constraint_type.upper()
    pocket_cutoff = float(args.pocket_cutoff)
    ncore = args.ncore

    ## => Run the Glide Docking Job <= ##

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    gridfile = os.path.abspath(gridfile)
    ligandfile = os.path.abspath(ligandfile)
    refligand = os.path.abspath(refligand)

    start_dir = os.getcwd()
    work_dir = os.path.join('docking', outdir)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    options = {
        'GRIDFILE' : gridfile,
        'LIGANDFILE' : ligandfile,
        'PRECISION' : precision,
        'SAMPLE_N_INVERSIONS' : True,
        'SAMPLE_RINGS' : True,
        'AMIDE_MODE' : 'penal',
        'EPIK_PENALTIES' : True,
        'COMPRESS_POSES' : True,
        'KEEP_SUBJOB_POSES' : True
    }

    if constraint_type == 'CORE':
        if refligand == None:
            raise Exception("The reference ligand file was not set for running CORE Constraint")

        options['CORE_RESTRAIN'] = True
        options['CORE_POS_MAX_RMSD'] = 0.10
        options['CORE_DEFINITION'] = "mcssmarts"
        options['CORE_FILTER'] = False
        options['USE_REF_LIGAND'] = True
        options['REF_LIGAND_FILE'] = refligand
    
    elif constraint_type == 'SHAPE':
        if refligand == None:
            raise Exception("The reference ligand file was not set for running SHAPE Constraint")

        options['SHAPE_RESTRAIN'] = True
        options['USE_REF_LIGAND'] = True
        options['SHAPE_REF_LIGAND_FILE'] = refligand
        options['REF_LIGAND_FILE'] = refligand

    elif constraint_type != 'NONE':
        raise Exception(f"The constraint {constraint_type} is not an available option!")

    if precision not in ['HTVS', 'SP', 'XP']:
        raise Exception(f"Precision type {precision} is not an available option!")

    dock_job = glide.Dock(options)

    dock_input = 'dockjob.inp'
    dock_job.writeSimplified(dock_input)

    subprocess.Popen([f'{schrodinger_path}/glide', dock_input, '-WAIT', '-HOST', f'localhost:{ncore}', '-OVERWRITE']).wait()

    pocketeer('dockjob_pv.maegz', 'pocket_pv.maegz', pocket_cutoff)

    os.chdir(start_dir)