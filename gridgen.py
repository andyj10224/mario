import subprocess, sys, os, argparse
import schrodinger.structure as structure
import schrodinger.structutils.analyze as analyze
import schrodinger.protein.findhets as findhets
import schrodinger.application.glide.glide as glide

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger Glide Gridgen Job on prepared ligands')

    parser.add_argument('protein_file', help='[string] The path to the prepared protein file')
    parser.add_argument('pdbid', help='[string] The pdbid of the protein being prepared')

    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])

    protein_file = args.protein_file
    pdbid = args.pdbid
    ncore = args.ncore

    ## ==> Run the Gridgen Docking Job <== ##

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Make the directory to hold the reference ligands
    ref_dir = os.path.dirname(protein_file)
    # Base name of the protein
    proteinbase = os.path.splitext(os.path.basename(protein_file))[0]

    protein_system = structure.StructureReader.read(protein_file)
    ligand_lists = findhets.find_hets(protein_system, include_metals=False, include_hydrogens=True)

    # Center of mass of all ligands in pdb file
    ligand_coms = []

    for n, atoms in enumerate(ligand_lists):
        ligand_st = protein_system.extract(atoms)
        ligand_obj = analyze.Ligand(ligand_st)
        ligand_coms.append(ligand_obj.centroid)
        ref_ligand_path = os.path.join(ref_dir, f'{proteinbase}_site_{n+1}_ligand.mae')
        structure.StructureWriter.write(ligand_st, ref_ligand_path)

    # Remove the ligands from the protein (and create the receptor file)
    atoms_to_remove = []

    for ligand in ligand_lists:
        atoms_to_remove.extend(ligand)

    protein_system.deleteAtoms(atoms_to_remove)
    protein_no_lig = os.path.join(ref_dir, f'{proteinbase}_no_ligand.mae')
    structure.StructureWriter.write(protein_system, protein_no_lig)

    protein_no_lig = os.path.abspath(protein_no_lig)

    options = {
        'FORCEFIELD' : 'OPLS_2005',
        'GLIDE_DIELCO' : 2.0,
        'GLIDECONS' : False,
        'INNERBOX' : [10, 10, 10],
        'OUTERBOX' : [30.0, 30.0, 30.0],
        'REC_MAECHARGES' : False,
        'RECEP_CCUT' : 0.25,
        'RECEP_VSCALE' : 1.0,
        'USECOMPMAE' : True,
        'RECEP_FILE' : protein_no_lig
    }

    start_dir = os.getcwd()
    work_dir = os.path.join('grids', pdbid)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    for n, com in enumerate(ligand_coms):
        options['JOBNAME'] = f'{proteinbase}_site_{n+1}'
        options['GRIDFILE'] = f'{proteinbase}_site_{n+1}.zip'
        options['GRID_CENTER'] = [com[0], com[1], com[2]]
        glide_job = glide.Gridgen(options)

        inp_file = f'{proteinbase}_site_{n+1}.inp'
        glide_job.writeSimplified(inp_file)

        subprocess.Popen([f'{schrodinger_path}/glide', inp_file, '-WAIT', '-HOST', f'localhost:{ncore}']).wait()
    
    os.chdir(start_dir)