import os
from subprocess import run, PIPE, Popen

from typing import Tuple

from .frame_counts import get_frames_count

import mdtraj

from .gmx_spells import merge_xtc_files

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

# Run a translation similar to gromacs trjconv '-trans' logic but spaced along the whole trajectory
# WARNING: This functions makes not sense alone but it may be useful combined with any further '-pbc -center'
def smooth_translation (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    translation : Tuple[float, float, float]
):
    # If the output trajectory file already exists we must stop here
    # The raw trjcat implementation will not join things to the end of it
    if os.path.exists(output_trajectory_filename):
        raise SystemExit('The output file already exists and overwrite is not supported for this command')
    # Split the translation among the number of frames in the trajectory
    frame_count = get_frames_count(input_structure_filename, input_trajectory_filename)
    frame_offsets = [ t / frame_count for t in translation ]
    # Now apply the translation on every frame
    # Rely on mdtraj for this step
    topology = mdtraj.load_topology(input_structure_filename)
    trajectory = mdtraj.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    frame_filename = '.frame.xtc'
    offsets = [0,0,0]
    print('\n')
    for f, frame in enumerate(trajectory):
        print(ERASE_PREVIOUS_LINE)
        print('Frame ' + str(f))
        offsets[0] += frame_offsets[0]
        offsets[1] += frame_offsets[1]
        offsets[2] += frame_offsets[2]
        xyz = frame.xyz[0]
        for atom in xyz:
            atom[0] += offsets[0]
            atom[1] += offsets[1]
            atom[2] += offsets[2]
        new_frame = mdtraj.Trajectory([xyz],
            topology=topology,
            time=frame.time,
            unitcell_lengths=frame.unitcell_lengths,
            unitcell_angles=frame.unitcell_angles
        )
        new_frame.save(frame_filename)
        merge_xtc_files(output_trajectory_filename, frame_filename)

    # Remove the residual files
    residual_files = [ frame_filename ]
    for residual_file in residual_files:
        os.remove(residual_file)