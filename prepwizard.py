import sys

import schrodinger.structure as structure
import schrodinger.application.prepwizard2.prepare as prepare
import schrodinger.application.prepwizard2.tasks as tasks

def prepare_pdb(pdbid : str) -> list:

    # Read in the pdb
    pdb_struct = prepare.retrieve_and_read_pdb(pdbid)
    
    ### => Preprocessing Input <= ###
    ppi = tasks.PreprocessInput()
    ppi.struct = pdb_struct

    # Preprocess Options
    # ppi.reference_structure = None
    ppi.preprocess_delete_far_waters = True
    ppi.preprocess_watdist = 5.0
    ppi.treat_metals = True
    ppi.assign_bond_orders = True
    ppi.use_ccd = True
    ppi.add_hydrogens = False
    ppi.readd_hydrogens = True
    ppi.add_terminal_oxygens = False
    ppi.treat_disulfides = True
    ppi.treat_glycosylation = False
    ppi.treat_palmitoylation = False
    # ppi.annotate_antibodies = True
    # ppi.antibody_cdr_scheme = annotation.DEFAULT_ANTIBODY_SCHEME
    ppi.selenomethionines = False
    ppi.fillsidechains = True
    ppi.fillloops = True
    # ppi.custom_fasta_file = None
    ppi.cap_termini = True
    ppi.run_epik = True
    ppi.idealize_hydrogen_tf = True

    # Epik Generate States Options
    gss = tasks.GenerateStatesSettings()
    gss.epik_pH = 7.4
    gss.epik_pHt = 0.1 # 2.0
    gss.max_states = 1
    gss.process_detected_ligands = True
    gss.process_metals_and_ions = True
    gss.process_non_water_solvents = False
    gss.process_other_hets = True

    ppi.generate_states_settings = gss

    ### => Optimize H Bond Input <= ###
    ohbi = tasks.OptimizeHBondInput()
    ohbi.samplewater = True
    ohbi.xtal = False
    ohbi.simplified_pH = 'neutral'
    ohbi.use_propka = True
    ohbi.propka_pH = 7.4
    ohbi.label_pkas = False
    ohbi.force_list = []
    ohbi.minimize_adj_h = False
    # ohbi.protassign_number_sequential_cycles = None
    # ohbi.protassign_max_cluster_size = None
    ohbi.idealize_hydrogen_tf = True

    ### => Cleanup/Force-Field Options <= ###
    ci = tasks.CleanupInput()
    ci.run_impref = True
    ci.force_field = "OPLS_2005"
    ci.rmsd = 0.30
    ci.fixheavy = False
    # Delete Waters by Distance
    ci.delete_far_waters = True
    ci.watdist = 0.0
    # Hydrogen Bonding
    ci.delete_nonbridging_waters = True
    ci.delwater_hbond_cutoff = 3

    ### => Preprocessing Workflow Settings <= ###
    ppwfi = tasks.PPWorkflowInput()
    ppwfi.struct = pdb_struct
    ppwfi.preprocess = ppi
    ppwfi.hbond = ohbi
    ppwfi.cleanup = ci
    ppwfi.do_preprocess = True
    ppwfi.do_hbond = True
    ppwfi.do_cleanup = True

    ### => Run the Protein Preparation Workflow <= ###
    ppwt = tasks.PPWorkflowTask()
    ppwt.input = ppwfi
    ppwt.backendMain()
    
    output_structs = ppwt.output.structs

    prepared_files = []

    for n, st in enumerate(output_structs):
        fname = f'{pdbid}_prepared_struct_{n}.mae'
        st.write(fname)
        prepared_files.append(fname)

    return prepared_files

if __name__ == '__main__':
    prepare_pdb(sys.argv[1])