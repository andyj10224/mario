import sys

import schrodinger.structure as structure

def apnetprep(posefile : str) -> None:
    for n, st in enumerate(structure.StructureReader(posefile)):
        name = st.title
        xyz = st.getXYZ()
        natom = st.atom_total
        net_charge = st.formal_charge
        properties = st.property
        if (n >= 1): gscore = properties.get('r_i_glide_gscore')

        print(name, net_charge)
        if (n >= 1): print(gscore)

if __name__ == '__main__':
    apnetprep(sys.argv[1])