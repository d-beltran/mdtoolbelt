# Functions powered by VMD
# Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. 
# http://www.ks.uiuc.edu/Research/vmd/

import os

from subprocess import run, PIPE, Popen

from .single_frame_getter import get_first_frame

# Set the script filename with all commands to be passed to vmd
commands_filename = '.commands.vmd'

# List all the vmd supported trajectory formats
vmd_supported_trajectory_formats = {'mdcrd', 'crd', 'dcd', 'xtc', 'trr'}

# Given a vmd supported topology with no coordinates and a single frame file, generate a pdb file
def vmd_to_pdb (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_structure_filename : str):

    single_frame_filename = get_first_frame(input_trajectory_filename, vmd_supported_trajectory_formats)

    # Prepare a script for VMD to run. This is Tcl language
    with open(commands_filename, "w") as file:
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Write the current topology in 'pdb' format
        file.write('$all frame first\n')
        file.write('$all writepdb ' + output_structure_filename + '\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        single_frame_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    os.remove(commands_filename)
# Set function supported formats
vmd_to_pdb.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'psf', 'parm'},
            'input_trajectory_filename': vmd_supported_trajectory_formats
        },
        'outputs': {
            'output_structure_filename': {'pdb'}
        }
    }
]