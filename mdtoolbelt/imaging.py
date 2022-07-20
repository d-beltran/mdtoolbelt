import os
from subprocess import run, PIPE, Popen

import mdtraj

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Run the gromacs trjconv '-pbc nojump' logic only over a specific subset of atoms
# WARNING: This may result in periodic boundary artifacts which may not be "legal"
# WARNING: You better know what you are doing
def selective_no_jump (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    selection : 'Selection'
):

    # Convert the selection to a ndx file gromacs can read
    selection_name = 'selection'
    ndx_selection = selection.to_ndx(selection_name)
    ndx_filename = '.selection.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)   

    # Run gromacs filtering the atom selection
    nojump_trajectory_filename = '.nojump.xtc'
    # DANI: truco temporal
    if not os.path.exists(nojump_trajectory_filename):
        p = Popen([
            "echo",
            selection_name,
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_structure_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            nojump_trajectory_filename,
            '-n',
            ndx_filename,
            '-pbc',
            'nojump',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # If the output file does not exists at this point then it means Gromacs failed for some reason
    if not os.path.exists(nojump_trajectory_filename):
        raise SystemExit('Something went wrong with Gromacs')

    # Generate also a pdb file to be used as structure further
    nojump_structure_filename = '.nojump.pdb'
    # DANI: truco temporal
    if not os.path.exists(nojump_structure_filename):
        p = Popen([
            "echo",
            selection_name,
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_structure_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            nojump_structure_filename,
            '-n',
            ndx_filename,
            '-dump',
            '0',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # If the output file does not exists at this point then it means Gromacs failed for some reason
    if not os.path.exists(nojump_structure_filename):
        raise SystemExit('Something went wrong with Gromacs')

    # Now merge the atoms of the 'nojump' trajectory with the rest of atoms of the original trajectory
    # Rely on mdtraj for this step
    original_topology = mdtraj.load_topology(input_structure_filename)
    original_trajectory = mdtraj.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    nojump_trajectory = mdtraj.iterload(nojump_trajectory_filename, top=nojump_structure_filename, chunk=1)
    frame_filename = '.frame.xtc'
    print('\n')
    for f, original_frame in enumerate(original_trajectory):
        print(ERASE_PREVIOUS_LINE)
        print('Frame ' + str(f))
        nojump_frame = next(nojump_trajectory)
        original_xyz = original_frame.xyz[0]
        nojump_xyz = nojump_frame.xyz[0]
        for i, atom_index in enumerate(selection.atom_indices):
            original_xyz[atom_index] = nojump_xyz[i]
        new_frame = mdtraj.Trajectory([original_xyz], topology=original_topology)
        new_frame.time = original_frame.time
        new_frame.save(frame_filename)
        merge_xtc_files(output_trajectory_filename, frame_filename)

    # Remove the residual files
    residual_files = [ ndx_filename, nojump_trajectory_filename, nojump_structure_filename, frame_filename ]
    for residual_file in residual_files:
        os.remove(residual_file)

# Join xtc files
def merge_xtc_files (current_file : str, new_file : str):
    # If the current file does nt exist then set the new file as the current file
    if not os.path.exists(current_file):
        os.rename(new_file, current_file)
        return
    # Run trjcat
    logs = run([
        "gmx",
        "trjcat",
        "-f",
        new_file,
        current_file,
        '-o',
        current_file,
        '-quiet'
    ],
    stdout=PIPE,
    stderr=PIPE
    ).stdout.decode()