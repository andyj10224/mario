import subprocess
import sys
import os
import time

import schrodinger.structure as structure
import schrodinger.structutils.analyze as analyze
import schrodinger.protein.findhets as findhets
import schrodinger.application.glide.glide as glide

keywords_list = """
ASLSTRINGS                = string_list(default=list()) # ASL strings for receptor scaling
CORE_FILTER               = boolean(default=False) # skip ligands that do not contain the core
CV_CUTOFF                 = float(default=0.0) # Coulomb-van der Waals energy cutoff used for final filtering
DIELMOD                   = option('rdiel', 'cdiel', default='rdiel') # type of dielectric to use: distance-dependent (rdiel) or constant (cdiel)
FLEXASL                   = string(default=None) # ASL defining flexible receptor atoms
FORCEFIELD                = string(default='OPLS_2005') # force field
GLIDE_CONS_FEAT_FILE      = string(default=None) # feature file name for constraints jobs
GLIDE_CONS_RMETCOORD      = float_list(default=list()) # sphere radii of Glide metal_coordination constraints
GLIDE_CONS_RNOEMAX        = float_list(default=list()) # maximum distances for Glide NOE constraints
GLIDE_CONS_RNOEMIN        = float_list(default=list()) # minimum distances for Glide NOE constraints
GLIDE_CONS_RPOS           = float_list(default=list()) # sphere radii of Glide positional constraints
GLIDE_CONS_XMETCOORD      = float_list(default=list()) # X-coordinates of Glide metal-coordination constraints
GLIDE_CONS_XNOE           = float_list(default=list()) # X-coordinates of targets for Glide NOE constraints
GLIDE_CONS_XPOS           = float_list(default=list()) # X-coordinates of Glide positional constraints
GLIDE_CONS_YMETCOORD      = float_list(default=list()) # Y-coordinates of Glide metal-coordination constraints
GLIDE_CONS_YNOE           = float_list(default=list()) # Y-coordinates of targets for Glide NOE constraints
GLIDE_CONS_YPOS           = float_list(default=list()) # Y-coordinates of Glide positional constraints
GLIDE_CONS_ZMETCOORD      = float_list(default=list()) # Z-coordinates of Glide metal-coordination constraints
GLIDE_CONS_ZNOE           = float_list(default=list()) # Z-coordinates of targets for Glide NOE constraints
GLIDE_CONS_ZPOS           = float_list(default=list()) # Z-coordinates of Glide positional constraints
GLIDE_DIELCO              = float(default=2.0, min=0.0, max=9999.9) # dielectric constant
GLIDE_ELEMENTS            = boolean(default=False) # run in "Glide Elements" mode
GLIDE_NTOTALCONS          = integer(default=0, min=0, max=10) # number of receptor atoms having constraints
GLIDE_NUMEXVOL            = integer(default=0, min=0) # number of receptor excluded-volume regions
GLIDE_NUMMETCOORDCONS     = integer(default=0, min=0) # number of receptor metal-coordination constraints
GLIDE_NUMMETCOORDSITES    = int_list(default=list()) # number of available coordination sites per metal-coordination constraint
GLIDE_NUMNOECONS          = integer(default=0, min=0) # number of receptor NOE constraints
GLIDE_NUMPOSITCONS        = integer(default=0, min=0) # number of receptor positional constraints
GLIDE_NUMUSEXVOL          = integer(default=0, min=0) # number of excluded-volume regions to use
GLIDE_RECEP_ASLSCALE      = boolean(default=False) # for Glide gridgen jobs, take per-atom charge and radius scaling factors from ASL specification
GLIDE_RECEP_MAESCALE      = boolean(default=False) # for Glide gridgen jobs, take per-atom charge and radius scaling factors from properties written by Maestro to the receptor file
GLIDE_REXVOL              = float_list(default=list()) # sphere radii of Glide excluded volumes
GLIDE_REXVOLIN            = float_list(default=list()) # inner sphere (max penalty) radii of Glide excluded volumes
GLIDE_XEXVOL              = float_list(default=list()) # X-coordinates of centers of Glide excluded volumes
GLIDE_YEXVOL              = float_list(default=list()) # Y-coordinates of centers of Glide excluded volumes
GLIDE_ZEXVOL              = float_list(default=list()) # Z-coordinates of centers of Glide excluded volumes
GLIDECONS                 = boolean(default=False) # use constraints
GLIDECONSATOMS            = int_list(default=list()) # receptor constraint atom list (H-bond and metal)
GLIDECONSNAMES            = string_list(default=list()) # constraint label list
GLIDECONSUSESYMATOMS      = bool_list(default=list()) # include symmetry-related receptor atoms for current H-bond constraint (default = all true)
GLIDERECEPTORSCALECHARGES = float_list(default=list()) # per-atom scale factors for receptor charges
GLIDERECEPTORSCALERADII   = float_list(default=list()) # per-atom scale factors for receptor radii
GLIDEXVOLNAMES            = string_list(default=list()) # excluded-volume label list
GRID_CENTER               = float_list(default=list(0.0, 0.0, 0.0)) # coordinates of the grid center
GRID_CENTER_ASL           = string(default=None) # ASL expression defining the grid center
GRIDFILE                  = string(default=None) # path to grid (.grd or .zip) file
HBOND_ACCEP_HALO          = boolean(default=False) # include halogens as possible H-bond acceptors in scoring
HBOND_CONSTRAINTS         = string_list(default=list()) # HBOND_CONSTRAINTS "hbond1 47", "hbond2 55"
HBOND_DONOR_AROMH         = boolean(default=False) # include aromatic H as a possible H-bond donor in scoring
HBOND_DONOR_AROMH_CHARGE  = float(default=0.0) # count aromatic H as a donor if its partial charge exceeds this value
HBOND_DONOR_HALO          = boolean(default=False) # include halogens as possible H-bond donors in scoring
INNERBOX                  = int_list(default=list(10, 10, 10)) # size of the box bounding the possible placements for the ligand centroid
JOBNAME                   = string(default='impact') # job name used for job control and as a filename prefix
LIGAND_INDEX              = integer(default=None, min=1) # index of the ligand entry to use for computing default GRID_CENTER and OUTERBOX
LIGAND_MOLECULE           = integer(default=None, min=1) # index of the ligand molecule to be removed
METAL_CONSTRAINTS         = string_list(default=list()) # METAL_CONSTRAINTS "metal1 64", ...
METCOORD_CONSTRAINTS      = string_list(default=list()) # METCOORD_CONSTRAINTS "label <ncoordsites>", ...
METCOORD_SITES            = string_list(default=list()) # METCOORD_SITES "<x1> <y1> <z1> <r1>", "<x2> <y2> <z2> <r2>", ...; Note that the list of "<x> <y> <z> <r>" tuples in a single METCOORD_SITES specification covers *all* of the specified METCOORD_CONSTRAINTS.  Thus the number of such tuples is equal to the *sum* of all the <ncoordsites> specifications.
NOE_CONSTRAINTS           = string_list(default=list()) # NOE_CONSTRAINTS "<label> <x> <y> <z> <rmin> <rmax>", ...
OUTERBOX                  = float_list(default=list(30.0, 30.0, 30.0)) # outer (grid) box size in Angstroms
OUTPUTDIR                 = string(default=None) # if present in gridgen, overrides directory path from GRIDFILE
PAIRDISTANCES             = float_list(default=None) # user-selected bond constraint distances
PEPTIDE                   = boolean(default=False) # use grid and sampling settings optimized for polypeptides
POSIT_CONSTRAINTS         = string_list(default=list()) # POSIT_CONSTRAINTS "<label> <x> <y> <z> <r>", ...
REC_MAECHARGES            = boolean(default=False) # use charges from receptor Maestro file instead of those from the force field
RECEP_CCUT                = float(default=0.25, min=0.0) # charge cutoff used to determine whether to use vdW scaling of receptor atoms
RECEP_FILE                = string(default=None) # receptor file name
RECEP_VSCALE              = float(default=1.0, min=0.0) # vdW scaling for receptor atoms (see also RECEP_CCUT)
REF_LIGAND_FILE           = string(default=None) # Glide reference ligand file name
SUBSTRATE_PENAL_FILE      = string(default='') # File listing the grid-cell coordinates and penalty values for substrate-envelope jobs
USECOMPMAE                = boolean(default=False) # write compressed output Maestro file (defaults to true for Glide jobs)
USEFLEXASL                = boolean(default=False) # specify movable atoms by ASL in Glide input file
USEFLEXMAE                = boolean(default=False) # specify movable atoms by property in input .mae file
"""

def make_grid(protein_file : str, pdbid : str) -> list:
    """
    The function used to call the pipeline, the central information flow hub.

    Args:
        protein_file (str) : The filepath of the prepared protein, located in prepwizard directory
        pdbid (str) : The pdbid of the protein we are working with (to help with naming and indexing)

    Returns:
        grid_files (list[str]) : The filepath of the protein grid files
    """

    schrodinger_path = os.environ.get('SCHRODINGER')
    if schrodinger_path is None: raise Exception("Environment variable $SCHRODINGER is not set.")

    proteinbase = os.path.splitext(os.path.split(protein_file)[-1])[0]

    protein_system = structure.StructureReader.read(protein_file)
    ligand_lists = findhets.find_hets(protein_system, include_metals=False, include_hydrogens=True)

    # Center of mass of all ligands in pdb file
    ligand_coms = []

    for n, atoms in enumerate(ligand_lists):
        ligand_st = protein_system.extract(atoms)
        ligand_obj = analyze.Ligand(ligand_st)
        ligand_coms.append(ligand_obj.centroid)

    # Remove the ligands from the protein (and create the receptor file)
    atoms_to_remove = []

    for ligand in ligand_lists:
        atoms_to_remove.extend(ligand)

    protein_system.deleteAtoms(atoms_to_remove)
    protein_no_lig = f'{os.path.splitext(protein_file)[0]}_no_ligand.mae'
    structure.StructureWriter.write(protein_system, protein_no_lig)

    protein_no_lig = os.path.abspath(protein_no_lig)

    options = {}
    options['ASLSTRINGS'] = []
    options['CORE_FILTER'] = False
    options['CV_CUTOFF'] = 0.0
    options['DIELMOD'] = 'rdiel'
    options['FLEXASL'] = None
    options['FORCEFIELD'] = 'OPLS_2005'
    options['GLIDE_CONS_FEAT_FILE'] = None
    options['GLIDE_CONS_RMETCOORD'] = []
    options['GLIDE_CONS_RNOEMAX'] = []
    options['GLIDE_CONS_RNOEMIN'] = []
    options['GLIDE_CONS_RPOS'] = []
    options['GLIDE_CONS_XMETCOORD'] = []
    options['GLIDE_CONS_XNOE'] = []
    options['GLIDE_CONS_XPOS'] = []
    options['GLIDE_CONS_YMETCOORD'] = []
    options['GLIDE_CONS_YNOE'] = []
    options['GLIDE_CONS_YPOS'] = []
    options['GLIDE_CONS_ZMETCOORD'] = []
    options['GLIDE_CONS_ZNOE'] = []
    options['GLIDE_CONS_ZPOS'] = []
    options['GLIDE_DIELCO'] = 2.0
    options['GLIDE_ELEMENTS'] = False
    options['GLIDE_NTOTALCONS'] = 0
    options['GLIDE_NUMEXVOL'] = 0
    options['GLIDE_NUMMETCOORDCONS'] = 0
    options['GLIDE_NUMMETCOORDSITES'] = []
    options['GLIDE_NUMNOECONS'] = 0
    options['GLIDE_NUMPOSITCONS'] = 0
    options['GLIDE_NUMUSEXVOL'] = 0
    options['GLIDE_RECEP_ASLSCALE'] = False
    options['GLIDE_RECEP_MAESCALE'] = False
    options['GLIDE_REXVOL'] = []
    options['GLIDE_REXVOLIN'] = []
    options['GLIDE_XEXVOL'] = []
    options['GLIDE_YEXVOL'] = []
    options['GLIDE_ZEXVOL'] = []
    options['GLIDECONS'] = False
    options['GLIDECONSATOMS'] = []
    options['GLIDECONSNAMES'] = []

    options['GLIDEXVOLNAMES'] = []
    
    options['HBOND_ACCEP_HALO'] = False
    options['HBOND_CONSTRAINTS'] = []
    options['HBOND_DONOR_AROMH'] = False
    options['HBOND_DONOR_AROMH_CHARGE'] = 0.0
    options['HBOND_DONOR_HALO'] = False
    options['INNERBOX'] = [10, 10, 10]

    options['METAL_CONSTRAINTS'] = []
    options['METCOORD_CONSTRAINTS'] = []
    options['METCOORD_SITES'] = []
    options['NOE_CONSTRAINTS'] = []
    options['OUTERBOX'] = [30.0, 30.0, 30.0]
    options['PEPTIDE'] = False
    options['POSIT_CONSTRAINTS'] = []
    options['REC_MAECHARGES'] = False
    options['RECEP_CCUT'] = 0.25

    options['RECEP_VSCALE'] = 1.0
    options['USECOMPMAE'] = True
    options['USEFLEXASL'] = False
    options['USEFLEXMAE'] = False

    options['RECEP_FILE'] = protein_no_lig

    grid_files = []

    start_dir = os.getcwd()
    work_dir = os.path.join('grids', pdbid)
    if not os.path.isdir(work_dir): os.makedirs(work_dir)
    os.chdir(work_dir)

    for n, com in enumerate(ligand_coms):
        options['JOBNAME'] = f'{proteinbase}_site_{n+1}_grid'
        options['GRIDFILE'] = f'{proteinbase}_site_{n+1}_grid.zip'
        options['GRID_CENTER'] = [com[0], com[1], com[2]]
        glide_job = glide.Gridgen(options)

        inp_file = f'{proteinbase}_site_{n+1}_grid.inp'
        glide_job.writeSimplified(inp_file)

        grid_files.append(os.path.join(work_dir, options['GRIDFILE']))

        if os.path.isfile(options['GRIDFILE']):
            raise Exception(f"{options['GRIDFILE']} already exists from a previous Glide job. \
                                Please remove the file to continue.")

        os.system(f'{schrodinger_path}/glide {inp_file}')

        while not os.path.isfile(options['GRIDFILE']): time.sleep(1.0)
    
    os.chdir(start_dir)
    return grid_files

if __name__ == '__main__':
    make_grid(sys.argv[1], sys.argv[2])