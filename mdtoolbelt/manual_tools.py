# Subset of tools used alone and directly from the console

import sys
import os
from subprocess import run, PIPE

# This tool allows you to set the chain of all atoms in a selection
# This is powered by VMD and thus the selection lenguage must be the VMD's
# Arguments are as follows:
# 1 - Input pdb filename
# 2 - Atom selection
# 3 - Chain letter
# 4 - Output pdb filename (Optional. Input filename by default)
# e.g. python chainer.py example.pdb 'resname POPS' 'M'
def chainer (
    input_pdb_filename : str,
    atom_selection : str,
    chain_letter : str,
    output_pdb_filename : str = None
):

    if not output_pdb_filename:
        output_pdb_filename = input_pdb_filename

    # Check the file exists
    if not os.path.exists(input_pdb_filename):
        raise SystemExit('ERROR: The file does not exist')

    # Set he path to a script with all commands needed for vmd to parse the topology file
    commands_filename = '.commands.vmd'

    # Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
    # In addition, if chains are missing, this script asigns chains by fragment
    # Fragments are atom groups which are not connected by any bond
    with open(commands_filename, "w") as file:
        # Select the specified atoms and set the specified chain
        file.write('set atoms [atomselect top "' + atom_selection + '"]\n')
        file.write('$atoms set chain ' + chain_letter + '\n')
        # Write the current topology in 'pdb' format
        file.write('set all [atomselect top "all"]\n')
        file.write('$all frame first\n')
        file.write('$all writepdb ' + output_pdb_filename + '\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_pdb_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    # Remove the vmd commands file
    os.remove(commands_filename)