# mario

Mario is a Schrodinger software wrapper that allows one to easily get geometries of protein-ligand docking poses
given a pdb code and a ligand geometry file (sdf).

## Using the code

Having the Schrodinger Software is a prerequisite. Call this code
`$SCHRODINGER/run main.py pdbid ligandfile.sdf`

IMPORTANT: If you do not have internet access, make sure that `pdbid.pdb`, with the exact name and format, is in the directory.

For example, `$SCHRODINGER/run main.py 4mxo S985ligand.sdf`