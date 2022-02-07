# Subset of tools used alone and directly from the console

import sys
import os
from subprocess import run, PIPE

from typing import Optional

# This tool allows you to set the chain of all atoms in a selection
# This is powered by VMD and thus the selection lenguage must be the VMD's
# Arguments are as follows:
# 1 - Input pdb filename
# 2 - Atom selection (All atoms by defualt)
# 3 - Chain letter (May be the flag 'fragment', which is the default indeed)
# 4 - Output pdb filename (Input filename by default)
def chainer (
    input_pdb_filename : str,
    atom_selection : Optional[str] = None,
    chain_letter : Optional[str] = None,
    output_pdb_filename : Optional[str] = None
):

    # If no atom selection is provided then all atoms are selected
    if not atom_selection:
        atom_selection = 'all'

    # If no chain letter is provided then the flag 'fragment' is used
    if not chain_letter:
        chain_letter = 'fragment'

    # If no output filename is provided then use input filename as output filename
    if not output_pdb_filename:
        output_pdb_filename = input_pdb_filename

    # Check the file exists
    if not os.path.exists(input_pdb_filename):
        raise SystemExit('ERROR: The file does not exist')

    # Set he path to a script with all commands needed for VMD to parse the topology file
    # The scripts is in Tcl lenguage
    commands_filename = '.commands.vmd'
       
    with open(commands_filename, "w") as file:
        # Select the specified atoms and set the specified chain
        file.write('set atoms [atomselect top "' + atom_selection + '"]\n')
        # In case chain letter is not a letter but the 'fragment' flag, asign chains by fragment
        # Fragments are atom groups which are not connected by any bond
        if chain_letter == 'fragment':
            # Get all different chain names
            file.write('set chains_sample [lsort -unique [${atoms} get chain]]\n')
            # Set letters in alphabetic order
            file.write('set letters "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"\n')
            # Get the number of fragments
            file.write('set fragment_number [llength [lsort -unique -integer [${atoms} get fragment]]]\n')
            # For each fragment, set the chain of all atoms which belong to this fragment alphabetically
            # e.g. fragment 0 -> chain A, fragment 1 -> chain B, ...
            file.write('for {set i 0} {$i <= $fragment_number} {incr i} {\n')
            file.write('	set fragment_atoms [atomselect top "fragment $i"]\n')
            file.write('	$fragment_atoms set chain [lindex $letters $i]\n')
            file.write('}\n')
            # Otherwise, set the specified chain
        else:
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