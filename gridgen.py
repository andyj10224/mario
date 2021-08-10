import subprocess, sys, os, argparse
import schrodinger.structure as structure
import schrodinger.structutils.analyze as analyze
import schrodinger.protein.findhets as findhets
import schrodinger.application.glide.glide as glide

if __name__ == '__main__':

    ## ==> Read in the arguments <== ##

    parser = argparse.ArgumentParser(description='Runs a Schrodinger Glide Gridgen Job on prepared ligands')

    parser.add_argument('pdbid', help='[string] The pdbid of the protein being prepared')
    parser.add_argument('ligid', help='[string] The id of the ligand at the binding site of the protein')

    parser.add_argument('--ncore', help='[int] Number of CPU cores to run the job on', default=1)

    args = parser.parse_args(sys.argv[1:])
    pdbid = args.pdbid
    ligid = args.ligid.upper()
    ncore = args.ncore

    ## ==> Run the Gridgen Docking Job <== ##

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    # Make sure the prepwizard job has already been run
    protein_file = f'prepwizard/{pdbid}/{pdbid}.mae'
    if not os.path.isfile(protein_file):
        raise Exception(f"The file prepwizard/{pdbid}/{pdbid}.mae has not been found. Please prepare your protein before calling Gridgen.")

    # The protein gets converted into a Schrodinger file, and which atoms represent ligands are determined
    protein_system = structure.StructureReader.read(protein_file)
    ligand_lists = findhets.find_hets(protein_system, include_metals=False, include_hydrogens=True)

    # Center of mass of the docking ligand in the pdb file
    ligand_com = [0.0, 0.0, 0.0]

    # Find the docking site ligand (match with ligid)
    docking_site_atoms = []
    for atoms in ligand_lists:
        list_ligid = protein_system.atom[atoms[0]].pdbres.strip()
        if list_ligid != ligid: continue
        ligand_st = protein_system.extract(atoms)
        ligand_obj = analyze.Ligand(ligand_st)
        ligand_com = ligand_obj.centroid
        ref_ligand_path = os.path.join('prepwizard', pdbid, f'{ligid}.mae')
        structure.StructureWriter.write(ligand_st, ref_ligand_path)
        docking_site_atoms = atoms
        break

    if len(docking_site_atoms) == 0:
        raise Exception(f"Docking site ligand, with the id {ligid.upper()} cannot be found in the pdb file!!!")
    
    # Remove the ligand at the binding site, and write the structure to a file
    protein_system.deleteAtoms(docking_site_atoms)
    protein_clear_docking_site = os.path.join('prepwizard', pdbid, f'{pdbid}_clear_docking_site.mae')
    structure.StructureWriter.write(protein_system, protein_clear_docking_site)

    protein_no_lig = os.path.abspath(protein_clear_docking_site)

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
        'RECEP_FILE' : protein_no_lig,
        'JOBNAME' : pdbid,
        'GRIDFILE' : f'{pdbid}.zip',
        'GRID_CENTER' : [ligand_com[0], ligand_com[1], ligand_com[2]]
    }

    # Change to grids directory
    start_dir = os.getcwd()
    work_dir = os.path.join('grids', pdbid)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    # Write the input file
    glide_job = glide.Gridgen(options)
    inp_file = f'{pdbid}.inp'
    glide_job.writeSimplified(inp_file)

    # Run the gridgen job
    gridgen_job = subprocess.Popen([f'{schrodinger_path}/glide', inp_file, '-WAIT', '-HOST', f'localhost:{ncore}', '-OVERWRITE'])
    gridgen_job.wait()

    # Move back to starting directory
    os.chdir(start_dir)

    if gridgen_job.returncode != 0:
        raise Exception("The grid generation job failed.")