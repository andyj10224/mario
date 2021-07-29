# MARIO

MARIO is a Drug Discovery Software Pipeline that integrates many different types of software, allowing the automation of the data transfer pathway from raw data to free energy calculations.

## Using the code

Having the Schrodinger Software is a prerequisite. The path to the Schrodinger suite must be set as an environmental variable `export SCHRODINGER="/path/to/schrodinger/suite"`

Call this code `python driver.py pdbid ligandfile` to run the pipeline.

IMPORTANT: Make sure that the pdb file is found in the `pdbs` directory as `pdbs/{pdbid}.pdb`, and the ligand file is found in the `ligands` directory as `ligands/{ligandfile}.sdf`.

Example call (provided in the directory), `python driver.py 4mxo S985`

## Other options

`--ncore` => Sets the number of CPU cores for parallel performance (Default: 1)

`--retrieve_pdb` => Allows the code to retrieve the pdb file online if the file does not exist at `pdbs/{pdbid}.pdb`

`--docking_precision` => The level of precision to run the docking job [HTVS, SP, or XP] (Default: SP)

`--constraint_type` => What kind of constraint to use for the docking [NONE, CORE, or SHAPE] (Default: NONE)

`--run_mmgbsa` => Allows you to automatically run MMGBSA on all the posefiles generated by the pipeline

`--run_apnet` => Allows you to run AP-Net-dG deltaG predictions on the posefiles generated by the pipeline (NOT YET IMPLEMENTED)

## Additional Information

For more information on how to run the driver or any of the processes:

Driver => `python driver.py -h`
Prepwizard => `SCHRODINGER/run prepwizard.py -h`
Ligprep => `SCHRODINGER/run ligprep.py -h` or `python ligprep.py -h`
Gridgen => `SCHRODINGER/run gridgen.py -h`
Docking => `SCHRODINGER/run dock.py -h`
MMGBSA => `SCHRODINGER/run mmgbsa.py -h` or `python mmgbsa.py -h`