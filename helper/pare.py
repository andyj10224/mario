import copy
from collections import defaultdict

from schrodinger.structutils import minimize, measure

PROPATOM_MONOMER = 'i_user_monomer_index'
PROP_DELETED_BONDS = 's_user_deleted_bonds'
PROP_ADDED_BONDS = 's_user_added_bonds'
PROP_FRAGMENTATION = 's_user_protein_treatment'
PROPBOND_FRAGMENT = 'b_user_frag_bond'

AAS = ['GLY ', 'ALA ', 'VAL ', 'LEU ', 'ILE ', 'MET ', 'PRO ',
       'SER ', 'THR ', 'CYS ', 'CYT ', 'CYX ', 'ASN ', 'GLN ',
       'PHE ', 'TYR ', 'TRP ', 'HIS ', 'HID ', 'HIE ', 'HIP ',
       'ASP ', 'ASH ', 'GLU ', 'GLH ', 'LYS ', 'LYN ', 'ARG ', 'ARN ']

# Dictionary of bonds to split for fragmenting peptide units
aa_fragment_atoms_bb = defaultdict(list, {
    'PRO ' : [ [' CD ', ' N  '] ],
    'ACE ' : [ [' CH3', ' C  '] ],
    'NME ' : [ [' CA ', ' N  '] ],
    } )
for AA in AAS:
    aa_fragment_atoms_bb[AA].extend([ [' CA ', ' N  '], [' CA ', ' C  '] ])

# Dictionary of bonds to split for fragmenting sidechains only
aa_fragment_atoms_bbsc = copy.deepcopy(aa_fragment_atoms_bb)
for AA in AAS:
    aa_fragment_atoms_bbsc[AA].extend([ [' CA ', ' CB '] ])

# Dictionary of bonds to split for fragmenting sidechains by functional group
aa_fragment_atoms_funcgroup_additional = {
    'SER ' : [ [' CB ', ' OG '] ],
    'THR ' : [ [' CB ', ' OG1'] ],
    'CYS ' : [ [' CB ', ' SG '] ],
    'CYT ' : [ [' CB ', ' SG '] ],
    'CYX ' : [ [' CB ', ' SG '] ],
    'ASN ' : [ [' CB ', ' CG '] ],
    'ASP ' : [ [' CB ', ' CG '] ],
    'ASH ' : [ [' CB ', ' CG '] ],
    'GLN ' : [ [' CG ', ' CD '] ],
    'GLU ' : [ [' CG ', ' CD '] ],
    'GLH ' : [ [' CG ', ' CD '] ],
    'PHE ' : [ [' CB ', ' CG '] ],
    'TYR ' : [ [' CB ', ' CG '] , [' CZ ', ' OH '] ],
    'TRP ' : [ [' CB ', ' CG '] ],
    'HIS ' : [ [' CB ', ' CG '] ],
    'HIE ' : [ [' CB ', ' CG '] ],
    'HID ' : [ [' CB ', ' CG '] ],
    'HIP ' : [ [' CB ', ' CG '] ],
    'LYS ' : [ [' CE ', ' NZ '] ],
    'LYN ' : [ [' CE ', ' NZ '] ],
    'ARG ' : [ [' CD ', ' NE '] ],
    'ARN ' : [ [' CD ', ' NE '] ],
    }
aa_fragment_atoms_funcgroup = copy.deepcopy(aa_fragment_atoms_bb)
for AA in aa_fragment_atoms_funcgroup_additional:
    aa_fragment_atoms_funcgroup[AA].extend(aa_fragment_atoms_funcgroup_additional[AA])

# Dictionary of bonds to split for fragmenting by functional group
aa_fragment_atoms_allalkyl_additional = {
    'VAL ' : [ [' CB ', ' CG1'], [' CB ', ' CG2'] ],
    'LEU ' : [ [' CB ', ' CG '], [' CG ', ' CD1'], [' CG ', ' CD2'] ],
    'ILE ' : [ [' CB ', ' CG1'], [' CB ', ' CG2'], [' CG1', ' CD1'] ],
    'MET ' : [ [' CB ', ' CG '], [' CG ', ' SD '], [' SD ', ' CE '] ],
    'PRO ' : [ [' CB ', ' CG '], [' CG ', ' CD '] ],
    'THR ' : [ [' CB ', ' CG2'] ],
    'GLN ' : [ [' CB ', ' CG '] ],
    'GLU ' : [ [' CB ', ' CG '] ],
    'GLH ' : [ [' CB ', ' CG '] ],
    'LYS ' : [ [' CB ', ' CG '], [' CG ', ' CD '], [' CD ', ' CE '] ],
    'LYN ' : [ [' CB ', ' CG '], [' CG ', ' CD '], [' CD ', ' CE '] ],
    'ARG ' : [ [' CB ', ' CG '], [' CG ', ' CD '] ],
    'ARN ' : [ [' CB ', ' CG '], [' CG ', ' CD '] ],
    }
aa_fragment_atoms_allalkyl = copy.deepcopy(aa_fragment_atoms_funcgroup)
for AA in AAS:
    aa_fragment_atoms_allalkyl[AA].append([' CA ', ' CB '])
for AA in aa_fragment_atoms_allalkyl_additional:
    aa_fragment_atoms_allalkyl[AA].extend(aa_fragment_atoms_allalkyl_additional[AA])

# Dictionary of atoms that must be retained as capping groups
aa_capping_atoms_bbsc = defaultdict(list, {('PRO ', ' N  ') : [' CD ']})
for AA in AAS:
    aa_capping_atoms_bbsc[(AA, ' CB ')].append(' CA ')

aa_capping_atoms_allalkyl = defaultdict(list, {
    ('MET ', ' SD ') : [ ' CG ', ' CE ' ],
    ('PRO ', ' N  ') : [ ' CD ' ],
    ('SER ', ' OG ') : [ ' CB ' ],
    ('THR ', ' OG1') : [ ' CB ' ],
    ('CYS ', ' SG ') : [ ' CB ' ],
    ('CYT ', ' SG ') : [ ' CB ' ],
    ('CYX ', ' SG ') : [ ' CB ' ],
    ('ASN ', ' CG ') : [ ' CB ' ],
    ('ASP ', ' CG ') : [ ' CB ' ],
    ('ASH ', ' CG ') : [ ' CB ' ],
    ('GLN ', ' CD ') : [ ' CG ' ],
    ('GLU ', ' CD ') : [ ' CG ' ],
    ('GLH ', ' CD ') : [ ' CG ' ],
    ('PHE ', ' CG ') : [ ' CB ' ],
    ('TYR ', ' CG ') : [ ' CB ' ],
    ('TRP ', ' CG ') : [ ' CB ' ],
    ('HIS ', ' CG ') : [ ' CB ' ],
    ('HID ', ' CG ') : [ ' CB ' ],
    ('HIE ', ' CG ') : [ ' CB ' ],
    ('HIP ', ' CG ') : [ ' CB ' ],
    ('LYS ', ' NZ ') : [ ' CE ' ],
    ('LYN ', ' NZ ') : [ ' CE ' ],
    ('ARG ', ' NE ') : [ ' CD ' ],
    ('ARN ', ' NE ') : [ ' CD ' ],
    } )

aa_bb_atoms = [ ' N  ', ' C  ', ' O  ' ]

aa_polar_atoms = {
    'SER ' : [ ' OG ' ],
    'THR ' : [ ' OG1' ],
    'CYS ' : [ ' SG ' ],
    'CYT ' : [ ' SG ' ], # CYS with S-
    'CYX ' : [ ' SG ' ], # CYS in disulfide
    'TYR ' : [ ' OH ' ],
    'TYO ' : [ ' OH ' ], # TYR with O-
    'ASN ' : [ ' OD1', ' CG ', ' ND2' ],
    'ASP ' : [ ' OD1', ' CG ', ' OD2' ],
    'ASH ' : [ ' OD1', ' CG ', ' OD2' ], # neutral
    'GLN ' : [ ' OE1', ' CD ', ' NE2' ],
    'GLU ' : [ ' OE1', ' CD ', ' OE2' ], # neutral
    'GLH ' : [ ' OE1', ' CD ', ' OE2' ],
    'LYS ' : [ ' NZ ' ],
    'LYN ' : [ ' NZ ' ], # neutral
    'ARG ' : [ ' NE ' , ' CZ ', ' NH1', ' NH2' ],
    'ARN ' : [ ' NE ' , ' CZ ', ' NH1', ' NH2' ], # neutral
    }

aa_aromatic_atoms = {
    'PHE ' : [ ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ' ],
    'TYR ' : [ ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ' ],
    'TYO ' : [ ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ' ],
    'TRP ' : [ ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2' ],
    'HIS ' : [ ' CG ', ' ND1', ' CD2', ' CE2', ' NE2' ],
    'HID ' : [ ' CG ', ' ND1', ' CD2', ' CE2', ' NE2' ],
    'HIE ' : [ ' CG ', ' ND1', ' CD2', ' CE2', ' NE2' ],
    'HIP ' : [ ' CG ', ' ND1', ' CD2', ' CE2', ' NE2' ],
    }


def fragment(s, protein_treatment='functional_group', atom_pairs=None):
    """Return a copy of structure s with fragment information

    Return a copy of structure object s in which residue information has been modified to
    divide molecules into "fragments."  Each molecule is a different fragment, and within
    each molecule fragments are created by breaking bonds in any of three ways:   
    1) Bonds within s that have property 'b_user_frag_bond' set to True
    2) If protein_treatment is one of 'sidechain' 'functional_group', or 'all_alkyl', amino acid residues are
       fragmented as described below
    3) Bonds formed by pairs of atoms contained in the list atom_pairs

    Protein treatment options:
    'sidechain': all peptide units and sidechains are treated as separate fragments.
    'functional_group': same as 'sidechain', with the addition that side chains with aromatic and polar groups are
       fragmented further into polar, aromatic, and aliphatic fragments.
    'all_alkyl': same as 'functional_group', with the addition that all methylene units in side chains are treated as
       individual fragments.

    :type s: schrodinger.structure.structure.Structure
    :param s: Structure to be fragmented

    :type protein_treatment: string
    :param protein_treatment: If not the empty string, automatically fragment proteins as specified. Must be one of the
        the following: 'sidechain', 'functional_group', 'all_alkyl'.

    :type atom_pairs: list of length-2 list of integers or None
    :param atom_pairs: list of indices of pairs of atoms in bonds to be fragmented

    :rtype: schrodinger.structure.structure.Structure
    :return: fragmented structure
    """
    if atom_pairs is None:
        atom_pairs = []

    sc = s.copy()

    # Create original residue information properties
    for at in sc.atom:
        if 's_user_pdbres_1' in at.property:
            at.pdbres = at.property['s_user_pdbres_1']
        else:
            at.property['s_user_pdbres_1'] = at.pdbres
        if 'i_user_resnum_1' in at.property:
            at.resnum = at.property['i_user_resnum_1']
        else:
            at.property['i_user_resnum_1'] = at.resnum
        if 's_user_inscode_1' in at.property:
            at.inscode = at.property['s_user_inscode_1']
        else:
            at.property['s_user_inscode_1'] = at.inscode

    # Identify all residues
    residue_inscodes = defaultdict(list)
    for r in sc.residue:
        residue_inscodes[(r.chain, r.pdbres, r.resnum)].append(r.inscode)

    # Break bonds by atom pair
    for iat1, iat2 in atom_pairs:
        bo = sc.getBond(iat1, iat2) # Returns None if no bond
        if bo:
            bo.property[PROPBOND_FRAGMENT] = True

    if protein_treatment == 'sidechain' or protein_treatment == 'sc':
        aa_fragment_dict = aa_fragment_atoms_bbsc
    elif protein_treatment == 'functional_group':
        aa_fragment_dict = aa_fragment_atoms_funcgroup
    elif protein_treatment == 'all_alkyl':
        aa_fragment_dict = aa_fragment_atoms_allalkyl
    elif protein_treatment:
        protein_treatment = ''
    sc.property[PROP_FRAGMENTATION] = protein_treatment
    if protein_treatment:
        for r in sc.residue:
            if r.pdbres not in aa_fragment_dict:
                continue
            frag_atoms = aa_fragment_dict[r.pdbres]
            for at1p, at2p in frag_atoms:
                at1 = r.getAtomByPdbName(at1p)
                at2 = r.getAtomByPdbName(at2p)
                if at1 and at2:
                    bo = sc.getBond(at1, at2) # Returns None if no bond
                    if bo:
                        bo.property[PROPBOND_FRAGMENT] = True

    deleted_bonds = []        
    for bo in sc.bond:
        try:
            if bo.property[PROPBOND_FRAGMENT]:
                deleted_bonds.append((bo.atom1, bo.atom2, bo.order))
        except KeyError:
            pass
    for iat1, iat2, order in deleted_bonds:
        sc.deleteBond(iat1, iat2)

    # Identify residues in each fragment and set fragment type
    residues_to_fragments = defaultdict(list)
    for mol in sc.molecule:
        iat = 1
        while True:
            at1 = mol.atom[iat]
            if at1.element != 'H':
                break
            iat += 1
        if at1.getResidue().isStandardResidue():
            if at1.pdbname in aa_bb_atoms:
                fragtype = 'pep'
                for at in mol.atom:
                    if at.pdbname == ' N  ':
                        at1 = at
            elif protein_treatment == 'sidechain' or protein_treatment == 'sc':
                fragtype = 'sc'
            elif at1.pdbres in aa_polar_atoms and at1.pdbname in aa_polar_atoms[at1.pdbres]:
                fragtype = 'pol'
            elif at1.pdbres in aa_aromatic_atoms and at1.pdbname in aa_aromatic_atoms[at1.pdbres]:
                fragtype = 'aro'
            else:
                fragtype = 'alk'
        else:
            fragtype = 'lig'
        residues_to_fragments[(at1.chain, at1.pdbres, at1.resnum, at1.inscode)].append((mol.number, fragtype))

    # Create inscode identifiers for fragments that have non-unique residue info
    for rinfo in residues_to_fragments:
        if len(residues_to_fragments[rinfo]) > 1:
            addl_inscodes = []
            for imol, fragtype in residues_to_fragments[rinfo]:
                mol = sc.molecule[imol]
                chain, pdbres, resnum, inscode = rinfo
                for insc in 'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'.split():
                    if insc not in residue_inscodes[(chain, pdbres, resnum)] + addl_inscodes:
                        break
                if fragtype == 'lig':
                    fragname = f'{chain}_{pdbres.strip()}{resnum}{inscode.strip()}_{fragtype}{insc}'
                elif fragtype == 'pep':
                    if inscode == ' ':
                        fragname = f'{chain}_{resnum-1}-{resnum}_{fragtype}'
                    else:
                        fragname = f'{chain}_{resnum}-{resnum}{inscode}_{fragtype}'
                else:
                    fragname = f'{chain}_{pdbres.strip()}{resnum}{inscode.strip()}_{fragtype}'
                for at in mol.atom:
                    at.chain = chain
                    at.pdbres = pdbres
                    at.resnum = resnum
                    at.inscode = insc
                    at.property['s_user_pdbres_f'] = pdbres
                    at.property['i_user_resnum_f'] = resnum
                    at.property['s_user_inscode_f'] = insc
                    at.property['s_user_fragment_name'] = fragname 
                addl_inscodes.append(insc)
        else:
            imol, fragtype = residues_to_fragments[rinfo][0]
            mol = sc.molecule[imol]
            chain, pdbres, resnum, inscode = rinfo
            if fragtype == 'lig':
                fragname = f'{chain}_{pdbres.strip()}{resnum}{inscode.strip()}_{fragtype}'
            elif fragtype == 'pep':
                if inscode == ' ':
                    fragname = f'{chain}_{resnum-1}-{resnum}_{fragtype}'
                else:
                    fragname = f'{chain}_{resnum}-{resnum}{inscode}_{fragtype}'
            else:
                fragname = f'{chain}_{pdbres.strip()}{resnum}{inscode.strip()}_{fragtype}'
            for at in mol.atom:
                at.property['s_user_fragment_name'] = fragname

    sc.addBonds(deleted_bonds)
    for at1, at2, order in deleted_bonds:
        bo = sc.getBond(at1.index, at2.index) # Returns None if no bond
        if bo:
            bo.property[PROPBOND_FRAGMENT] = True

    return sc


def unfragment(s, restore_atom=False, restore_bond=False):
    """Return a copy of structure s with fragment information eliminated

    Returns a copy of structure object s with residue information that had been changed by function fragment
    restored to the original information

    :type s: schrodinger.structure.structure.Structure
    :param s: previously fragmented structure to be restored

    :type restore_atom: boolean
    :param restore_atom: If True, eliminate atom properties that were created by function fragment

    :type restore_bond: boolean
    :param restore_bond: If True, eliminate bond properties that were created by function fragment

    :rtype: schrodinger.structure.structure.Structure
    :return: restored structure
    """
    sc = s.copy()
    if PROP_FRAGMENTATION in sc.property:
        del sc.property[PROP_FRAGMENTATION]

    for at in sc.atom:
        if 's_user_pdbres_1' in at.property:
            at.pdbres = at.property['s_user_pdbres_1']
        if 'i_user_resnum_1' in at.property:
            at.resnum = at.property['i_user_resnum_1']
        if 's_user_inscode_1' in at.property:
            at.inscode = at.property['s_user_inscode_1']

    if restore_atom:
        for at in sc.atom:
            for prop in 's_user_pdbres_1 i_user_resnum_1 s_user_inscode_1 s_user_pdbres_f i_user_resnum_f s_user_inscode_f'.split():
                if prop in at.property:
                    del at.property[prop]

    if restore_bond:
        for bo in sc.bond:
            if PROPBOND_FRAGMENT in bo.property:
                del bo.property[PROPBOND_FRAGMENT]

    return sc


def addFullFragmentFromAtomList(s, atom_list=None, fa_set=None, h_dict=None, debug=False):
    """Add atom indices of new fragments

    Given a fragmented structure and a list of atom indices, identify the fragments containing one or
    more of the atom indices. All atoms in these fragments are added to pre-defined set fa_set, and
    h_dict is updated appropriately. If the addition of the fragment creates any uncapped fragment(s),
    add an appropriate capping fragment also. If addition of the fragment creates a situation where a
    single atom functions as the capping hydrogen for two or more fragments, add the fragment
    containing that atom also.

    :type s: schrodinger.structure.structure.Structure
    :param s: fragmented structure

    :type atom_list: list of integers or None
    :param atom_list: list of indices for atoms that determine the fragments to be added

    :type fa_set: set of integers or None
    :param fa_set: initial set of atom indices to which to add fragment atoms

    :type h_dict: dictionary (key:value type int:int) or None
    :param h_dict: if not None, dictionary where each key is an index of an atom which should be treated as a capping
        hydrogen and its value is the index of the fragment atom to which it is bonded
    
    :type debug: boolean
    :param debug: If True, print useful debugging information to stdout

    :rtype: tuple of length 2
    :return: updated set of fragment atom indices and h_dict
    """
    if atom_list is None:
        atom_list = []
    if fa_set is None:
        fa_set = set()
    if h_dict is None:
        h_dict = {}

    h_dict_new = copy.copy(h_dict)  # Copy h_dict to avoid changing to object directly
    f_set = set([s.atom[at].getResidue() for at in fa_set])  # Set of fragments (residues) will be converted to a set of atom indices at the end
    f_added = set([s.atom[at].getResidue() for at in atom_list])  # f_added: fragments added "recently"
    if debug:
        print(f'Adding residues to set:')
        for a1 in [at for f in f_added for at in f.atom]:
            print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')

    fragmentation = s.property[PROP_FRAGMENTATION]
    if fragmentation == 'sidechain' or fragmentation == 'sc':
        aa_capping_atoms = aa_capping_atoms_bbsc
    elif fragmentation == 'functional_group' or fragmentation == 'all_alkyl':
        aa_capping_atoms = aa_capping_atoms_allalkyl
    else:
        aa_capping_atoms = {}

    # Iterate through loop until no more fragments have been added; this needs to be iterative
    # because in some situations, adding some fragments may force the addition of other fragments
    while f_added:
        f_set |= f_added  # Add fragments to full set
        f_added_in_loop = []  # Fragments added in this iteration

        # Analyze all atoms that are bonded to atoms in the f_added fragments
        for r in f_added:
            for at in r.atom:
                for ba in at.bonded_atoms:
                    bar = ba.getResidue()
                    bai = ba.index
                    # If bonded atom already present in set: do nothing
                    if bar in f_set:
                        continue
                    # If fragment atom and bonded atom are pre-defined as a required cap:
                    #   add bonded fragment for processing in next iteration
                    elif ( ba.pdbname in aa_capping_atoms[(at.pdbres, at.pdbname)] or
                           ba.pdbname in [' CA ', ' CH3' ] and at.pdbname == ' C  ' or  # CH3 required for ACE residues
                           ba.pdbname == ' CA ' and at.pdbname == ' N  ' ):
                        f_added_in_loop.append(bar)
                        if debug:
                            print(f'- Extending to include capping fragment:')
                            for a1 in bar.atom:
                                print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')
                    # If bonded atom has been previously designated as a capping H atom for a different atom:
                    #   add bonded fragment for processing in next iteration
                    elif bai in h_dict_new and h_dict_new[bai] != at.index:
                        f_added_in_loop.append(bar)
                        del h_dict_new[bai]
                        if debug:
                            print(f'- Extending to include fragment from double-H-cap:')
                            for a1 in bar.atom:
                                print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')
                    # Otherwise: designate the bonded atom as a capping H atom
                    elif bai not in h_dict_new:
                        h_dict_new[bai] = at.index
                        if debug:
                            print(f'- Bonded atom added to H-cap list: {bai} {ba.pdbres} {ba.resnum} {ba.inscode} {ba.pdbname}')

        f_added = set(f_added_in_loop)

    # Convert added fragments (residues) to atom indices and add to fa_set
    fa_set_new = fa_set | set([at.index for r in f_set for at in r.atom])
    # Update h_dict_new by removing atoms from h_dict_new if they are now contained in fa_set
    for iat in fa_set_new:
        try:
            del h_dict_new[iat]
        except KeyError:
            pass
        
    return fa_set_new, h_dict_new


def removeFullFragmentFromAtomList(s, atom_list=None, fa_set=None, h_dict=None, debug=False):
    """Remove atom indices of fragments

    Given a fragmented structure and a list of atom indices, identify the fragments containing one or
    more of the atom indices. All atoms in these fragments are removed from pre-defined set fa_set, and
    h_dict is updated appropriately. If removal of the fragment creates any uncapped fragment(s), remove
    those fragments also. If removal of the fragment creates a situation where a single atom functions
    as the capping hydrogen for two or more fragments, remove those fragments also.

    :type s: schrodinger.structure.structure.Structure
    :param s: fragmented structure

    :type atom_list: list of integers or None
    :param atom_list: list of indices for atoms that determine the fragments to be removed

    :type fa_set: set of integers or None
    :param fa_set: initial set of atom indices from which to remove fragment atoms

    :type h_dict: dictionary (key:value type int:int) or None
    :param h_dict: dictionary where each key is an index of an atom which should be treated as a capping
        hydrogen and its value is the index of the fragment atom to which it is bonded
    
    :type debug: boolean
    :param debug: If True, print useful debugging information to stdout

    :rtype: tuple of length 2
    :return: updated set of fragment atom indices and h_dict
    """
    if atom_list is None:
        atom_list = []
    if fa_set is None:
        fa_set = set()
    if h_dict is None:
        h_dict = {}

    h_dict_new = copy.copy(h_dict)  # Copy h_dict to avoid changing to object directly
    f_set = set([s.atom[at].getResidue() for at in fa_set])  # Set of fragments (residues) will be converted to a set of atom indices at the end
    f_removed = set([s.atom[at].getResidue() for at in atom_list])  # f_removed: fragments removed "recently"
    if debug:
        print(f'Removing residues from set:')
        for a1 in [at for f in f_removed for at in f.atom]:
            print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')

    fragmentation = s.property[PROP_FRAGMENTATION]
    if fragmentation == 'sidechain' or fragmentation == 'sc':
        aa_capping_atoms = aa_capping_atoms_bbsc
    elif fragmentation == 'functional_group' or fragmentation == 'all_alkyl':
        aa_capping_atoms = aa_capping_atoms_allalkyl
    else:
        aa_capping_atoms = {}

    # Iterate through loop until no more fragments have been removed; this needs to be iterative
    # because in some situations, removing some fragments may force the removal of other fragments
    while f_removed:
        f_set -= f_removed  # Remove fragments from full set
        f_removed_in_loop = []  # Fragments removed in this iteration
        # Required to handle the undesirable situation where at atom ends up as a capping H for two different atoms
        # Note value type in h_dict_in_loop (list) is different from value type of h_dict_new (int)
        h_dict_in_loop = defaultdict(list)  

        # Analyze all atoms that are bonded to atoms in the f_removed fragments
        for r in f_removed:
            for at in r.atom:
                if debug:
                    print(f'Checking atom {at.index}')
                for ba in at.bonded_atoms:
                    bar = ba.getResidue()
                    bai = ba.index
                    if debug:
                        print(f'* bonded atom {bai}, {bar.pdbres} {bar.resnum} {bar.inscode} {ba.pdbname}')
                    # If bonded atom is in h_dict_new: remove it from h_dict_new
                    if bai in h_dict_new:
                        if debug:
                            print(f'  - Removing bonded atom {bai} from h_dict')
                        if h_dict_new[bai] != at.index:
                            print(f'Warning: Unusual situation detected: attempting to remove atom {bai} from h_dict due to removal of fragment containing bonded atom {at.index}, but h_dict value does not match: h_dict[{bai}] = {h_dict_new[bai]}')
                        del h_dict_new[bai]
                    # If bonded atom is already absent from set: do nothing
                    elif bar not in f_set:
                        if debug:
                            print(f'  - Doing nothing; not in set')
                        continue
                    # If fragment atom and bonded atom are pre-defined as a required cap:
                    #   add bonded fragment for processing in next iteration
                    elif ( at.pdbname in aa_capping_atoms[(ba.pdbres, ba.pdbname)] or
                           at.pdbname in [' CA ', ' CH3' ] and ba.pdbname == ' C  ' or  # CH3 required for ACE residues
                           at.pdbname == ' CA ' and ba.pdbname == ' N  ' ):
                        f_removed_in_loop.append(bar)
                        if debug:
                            print(f'  ** Will remove additional fragment due to cap removal:')
                            for a1 in bar.atom:
                                print(f'    {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')
                    # Otherwise: designate the fragment atom as a capping H atom for bonded atom
                    else:
                        h_dict_in_loop[at.index].append(bai)
                        if debug:
                            print(f'  - setting h_dict_in_loop: {at.index}; now at {h_dict_in_loop[at.index]}')

        # Check the atoms designated at H atoms to see if any was designated more than once
        # If so, add the bonded fragments for processing in next iteration
        # If not, update h_dict_new
        for iat in h_dict_in_loop:
            if len(h_dict_in_loop[iat]) > 1:
                bars = [s.atom[bai].getResidue() for bai in h_dict_in_loop[iat]]
                f_removed_in_loop.extend(bars)
                if debug:
                    print(f'  ** Removing additional fragment(s) to avoid double-H-cap:')
                    for bar in bars:
                        for a1 in bar.atom:
                            print(f'    {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')
            elif len(h_dict_in_loop[iat]) == 1:
                h_dict_new[iat] = h_dict_in_loop[iat][0]

        f_removed = set(f_removed_in_loop)

    return set([at.index for r in f_set for at in r.atom]), h_dict_new


def pareProtein(p, l=None, atom_list=None, nmax=0, debug=False):
    """Pare protein to fragments based on distance from the ligand

    Pare protein down to a smaller set of atoms, ensuring that only capped full fragments are present.
    All fragments containing any atom in atom_list are automatically retained in the pared set. If the
    pared set contains fewer than nmax atoms, fragments are added (prioritizing those nearest to the
    ligand) so that the number of atoms in the pared set does not exceed nmax.

    :type p: schrodinger.structure.structure.Structure
    :param p: fragmented structure representing the protein (may also include ligands). In p, fragments
        are defined as different residues. Use dimer.fragment before calling pareProtein to
        fragment it in a meaningful way.

    :type l: None or schrodinger.structure.structure.Structure
    :param l: if not None, structure representing the ligand (used to identify distances between protein fragments
        and the ligand)

    :type atom_list: list of integers or None
    :param atom_list: list of indices of the atoms determining fragments that will be automatically
        included in the pared set

    :type nmax: int
    :param nmax: maximum number of atoms in the pared set. This is a "soft" cutoff in the sense that
        it is not guaranteed that the number of atoms will not exceed nmax. Specifically, all of the
        atoms in all capped fragments obtained from parameter atom_list are retained, even if this
        number exceeds nmax.

    :type debug: boolean
    :param debug: If True, print useful debugging information to stdout

    :rtype: tuple of length 2
    :return: list of all atom indices in pared set, dictionary where keys are indices are atoms that
        should be treated as hydrogen caps
    """
    if atom_list is None:
        atom_list = []

    # Calculate minimum distances between fragments and ligand
    if nmax > 0 and l is not None:
        dmins = []
        for r in p.residue:
            dmin = 9999999
            for pat in r.atom:
                for lat in l.atom:
                    d = measure.measure_distance(pat, lat)
                    if d < dmin:
                        dmin = d
            dmins.append((dmin, r))
        dmins.sort(reverse = True)

    # Create initial fragmented atom set from supplied list
    fa_set, h_dict = addFullFragmentFromAtomList(p, atom_list = atom_list, debug = debug)

    # If atom set is still smaller than nmax atoms, build up atom set, closest fragment first
    at_count = len(fa_set) + len(h_dict)
    while at_count < nmax and dmins:
        dmin, r_close = dmins.pop()  # Next closest fragment
        if debug:
            print('* New fragment added to pared set:')
            for a1 in r_close.atom:
                print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')

        # Copy fa_set and h_dict in case we exceed nmax and need to revert to the last ones
        fa_set_copy = copy.copy(fa_set)
        h_dict_copy = copy.copy(h_dict)

        # Update atom set
        fa_set_copy, h_dict_copy = addFullFragmentFromAtomList(p, atom_list = [at.index for at in r_close.atom], fa_set = fa_set_copy, h_dict = h_dict_copy, debug = debug)
        at_count = len(fa_set_copy) + len(h_dict_copy)
        if debug:
            print(f'* Atom count: {at_count} ( {at_count - len(h_dict_copy)} / {len(h_dict_copy)} )')
            print(f'     fa_set: {["{} {} {}".format(r.pdbres, r.resnum, r.inscode) for r in set([p.atom[at].getResidue() for at in fa_set])]}\n     h_dict: {["{} {} {} {} {}".format(p.atom[iat].index, p.atom[iat].pdbres, p.atom[iat].resnum, p.atom[iat].inscode, p.atom[iat].pdbname) for iat in h_dict]}\n')

        if at_count <= nmax:
            fa_set = fa_set_copy
            h_dict = h_dict_copy

        else:  # Terminating condition: exceed max atom size
            if debug:
                print(f'* Max atoms exceeded: {at_count}; ignoring previous fragment update')

    if debug:
        print(f'* Final sets:')
        print(f'    fa_set: {fa_set}\n    h_dict: {h_dict}')
    return fa_set, h_dict


def pareProteinWithinDistance(p, l=None, atom_list=None, dmax=0.0, debug=False):
    """Pare protein to fragments based on distance from the ligand

    Pare protein down to a smaller set of atoms, ensuring that only capped full fragments are present.
    All fragments containing any atom in atom_list are automatically retained in the pared set. Any fragments that
    contain an atom within the cutoff distance dmax from any ligand atom are also included.

    :type p: schrodinger.structure.structure.Structure
    :param p: fragmented structure representing the protein (may also include ligands). In p, fragments
        are defined as different residues. Use dimer.fragment before calling pareProtein to
        fragment it in a meaningful way.

    :type l: None or schrodinger.structure.structure.Structure
    :param l: if not None, structure representing the ligand (used to identify distances between protein fragments
        and the ligand)

    :type atom_list: list of integers or None
    :param atom_list: list of indices of the atoms determining fragments that will be automatically
        included in the pared set

    :type dmax: float
    :param dmax: maximum distance for atoms in the pared set. A fragment will be retained if any of its atoms
        are within dmax Angstroms from any atom in the ligand.

    :type debug: boolean
    :param debug: If True, print useful debugging information to stdout

    :rtype: tuple of length 2
    :return: list of all atom indices in pared set, dictionary where keys are indices are atoms that
        should be treated as hydrogen caps
    """
    if atom_list is None:
        atom_list = []

    # Calculate minimum distances between fragments and ligand
    if dmax > 0 and l is not None:
        dmins = []
        for r in p.residue:
            dmin = 9999999
            for pat in r.atom:
                for lat in l.atom:
                    d = measure.measure_distance(pat, lat)
                    if d < dmin:
                        dmin = d
            dmins.append((dmin, r))
        dmins.sort(reverse = True)

    # Create initial fragmented atom set from supplied list
    fa_set, h_dict = addFullFragmentFromAtomList(p, atom_list = atom_list, debug = debug)

    # If atom set is still smaller than nmax atoms, build up atom set, closest fragment first
    dmin = 0
    while dmin <= dmax and dmins:
        dmin, r_close = dmins.pop()  # Next closest fragment
        if debug:
            print('* New fragment added to pared set:')
            for a1 in r_close.atom:
                print(f'   {a1.index} {a1.pdbres} {a1.resnum} {a1.inscode} {a1.pdbname}')

        # Update atom set
        fa_set, h_dict = addFullFragmentFromAtomList(p, atom_list = [at.index for at in r_close.atom], fa_set = fa_set, h_dict = h_dict, debug = debug)
        at_count = len(fa_set) + len(h_dict)
        if debug:
            print(f'* Atom count: {at_count} ( {at_count - len(h_dict)} / {len(h_dict)} )')
            print(f'     fa_set: {["{} {} {}".format(r.pdbres, r.resnum, r.inscode) for r in set([p.atom[at].getResidue() for at in fa_set])]}\n     h_dict: {["{} {} {} {} {}".format(p.atom[iat].index, p.atom[iat].pdbres, p.atom[iat].resnum, p.atom[iat].inscode, p.atom[iat].pdbname) for iat in h_dict]}\n')

    if debug:
        print(f'* Final sets:')
        print(f'    fa_set: {fa_set}\n    h_dict: {h_dict}')
    return fa_set, h_dict


def pareMinimize(s, atom_set=None, h_dict = None):
    """Pare and minimize system

    Reduce structure s to only the atoms in atom_set and h_dict. All atoms in h_dict are converted
    to H atoms. MacroModel is used to minimize only those H atoms using the default force field.

    :type s: schrodinger.structure.structure.Structure
    :param s: structure to be pared and minimized

    :type atom_set: set or None
    :param atom_set: set of atom indices to retain in s

    :type h_dict: dictionary or None
    :param h_dict: dictionary in which the keys are indices of atoms within s to be modified to H

    :rtype: tuple of length 2
    :return: pared, minimized structure and the renumbering dictionary to convert atom indices from the input s to the pared system
    """
    if atom_set is None:
        atom_set = set()
    if h_dict is None:
        h_dict = {}

    sc = s.copy()

    # Pare system to only atoms in atom_set and h_dict
    keep_atoms = atom_set | set(h_dict.keys())
    renumbering = sc.deleteAtoms([at.index for at in sc.atom if at.index not in keep_atoms], renumber_map = True)

    # Create new h_dict with proper renumbering
    h_dict_r = {}
    for iat in h_dict:
        h_dict_r[renumbering[iat]] = renumbering[h_dict[iat]]

    # Modify all atoms in h_dict
    # * change element to H
    # * change color to white
    # * remove any bonds other than to its non-H atom
    # * change PDB name; name is same as bonded atom, except 'C' replaced with 'H' and appending '0'
    # * retype
    for iat in h_dict_r:
        at = sc.atom[iat]  # Atom object of "H" atom
        at2 = sc.atom[h_dict_r[iat]]  # Atom object to which H is bonded

        at.element = 'H'
        at.color = 'white'
        bonded_atoms_to_delete = [bonded_atom for bonded_atom in at.bonded_atoms if bonded_atom != at2]
        for ba in bonded_atoms_to_delete:
            sc.deleteBond(at, ba)
        pdbn = at2.pdbname.replace('C', 'H', 1)
        pdbn = pdbn[1:] + pdbn[0]
        pdbn = pdbn.replace(' ', '0', 1)
        at.pdbname = pdbn[-1] + pdbn[:-1]
        at.retype()

    sc2 = sc.copy()

    # Minimize only newly created H atoms
    scmin = minimize.Minimizer(struct = sc)
    for iat in atom_set:
        scmin.addPosFrozen(renumbering[iat])
    scmin.minimize()

    # Copy coordinates only to new structure (to avoid creating extraneous properties from the minimization)
    sc2.setXYZ(sc.getXYZ())
    return sc2, renumbering