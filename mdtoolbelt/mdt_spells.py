import os
from subprocess import run, PIPE
from typing import List

import mdtraj as mdt
import numpy as np

from .formats import get_format
from .frame_counts import get_frames_count
from .gmx_spells import merge_xtc_files

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Set mdtraj supported formats
mdtraj_supported_structure_formats = {
    'pdb', 'pdb.gz' 'h5', 'lh5', 'prmtop', 'parm7', 'prm7', 'psf', 'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'
}
mdtraj_supported_trajectory_formats = {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos', 'mdcrd', 'xyz'}

# Use mdtraj 'mdconvert' command-line script (there is no python version for this tool apparently)
# Multiple files may be selected with bash syntax
def merge_and_convert_trajectories (
    input_trajectory_filenames : List[str],
    output_trajectory_filename : str
    ):

    # Run MDtraj
    logs = run([
        "mdconvert",
        "-o",
        output_trajectory_filename,
        *input_trajectory_filenames,
    ], stderr=PIPE).stderr.decode()
    # If output has not been generated then warn the user
    if not os.path.exists(output_trajectory_filename):
        print(logs)
        raise SystemExit('Something went wrong with MDTraj')
        
# NEVER FORGET: mdconvert does not handle mdcrd format
mdconvert_supported_formats = {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos'}
merge_and_convert_trajectories.format_sets = [
    {
        'inputs': {
            'input_trajectory_filenames': mdconvert_supported_formats
        },
        'outputs': {
            'output_trajectory_filename': mdconvert_supported_formats
        }
    },
]

# Merge and convert trajectories without the mdconvert command
# WARNING: This process is slow since we must iterate and merge each frame
# WARNING: This process is restricted to trr/xtc output files given the merger
def merge_and_convert_trajectories_alternative (
    input_structure_filename : str,
    input_trajectory_filenames : List[str],
    output_trajectory_filename : str
    ):

    # Assert we have input values
    if not input_structure_filename:
        raise SystemExit('ERROR: Missing input structure filenames')
    if not input_trajectory_filenames:
        raise SystemExit('ERROR: Missing input trajectory filenames')
    if not output_trajectory_filename:
        raise SystemExit('ERROR: Missing output trajectory filename')

    # Load the topology, which is used further
    topology = mdt.load_topology(input_structure_filename)

    # If the output trajectory filename matches any input file then we must rename the input filename to not overwrite it
    if output_trajectory_filename in input_trajectory_filenames:
        backup_filename = 'backup.' + output_trajectory_filename
        os.rename(output_trajectory_filename, backup_filename)
        repeated_input_filename_index = input_trajectory_filenames.index(output_trajectory_filename)
        input_trajectory_filenames[repeated_input_filename_index] = backup_filename

    # If the output trajectory file already exists at this point then we must stop here
    # The raw trjcat implementation will not join things to the end of it
    if os.path.exists(output_trajectory_filename):
        raise SystemExit('The output file already exists and overwrite is not supported for this functionality')
    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()
    # Iterate over the different input trajectory filenames
    frame_filename = '.frame.xtc'
    for input_trajectory_filename in input_trajectory_filenames:
        # Load the trajectory frame by frame
        trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
        for f, frame in enumerate(trajectory):
            # Update the current frame log
            print(ERASE_PREVIOUS_LINE)
            print('Frame ' + str(f))
            frame.save(frame_filename)
            # Join current frame to the output trajectory
            merge_xtc_files(output_trajectory_filename, frame_filename)
    # Remove the residual file
    # WARNING: It may not exist if the trajectory has 1 frame
    if os.path.exists(frame_filename):
        os.remove(frame_filename)
merge_and_convert_trajectories_alternative.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filenames': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    },
]

# Merge and convert trajectories without the mdconvert command
# WARNING: Note that this process is not memory efficient so beware the size of trajectories to be converted
# DEPRECATED: This function was meant to convert trajectories to mdcrd, which is not supported by mdconvert
# DEPRECATED: However the export to mdcrd in mdtraj does not allow to remove the periodic box, which may be a problem
# DEPRECTAED: Use VMD instead
def merge_and_convert_trajectories_unefficient (
    input_structure_filename : str,
    input_trajectory_filenames : List[str],
    output_trajectory_filename : str,
):
    print('WARNING: You are using a not memory efficient tool. If the trajectory is too big your system may not hold it.')
    # Load all trajectories together
    sample_trajectory = input_trajectory_filenames[0]
    if input_structure_filename:
        trajectory = mdt.load(sample_trajectory, top=input_structure_filename)
    else:
        trajectory = mdt.load(sample_trajectory)
    for extra_trajectory in input_trajectory_filenames[1:]:
        if input_structure_filename:
            extra = mdt.load(sample_trajectory, top=input_structure_filename)
        else:
            extra = mdt.load(sample_trajectory)
        trajectory = mdt.join([trajectory, extra], check_topology=False)
        
    # Write the whole trajectory
    trajectory.save(output_trajectory_filename)
    
merge_and_convert_trajectories_unefficient.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filenames': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_trajectory_formats
        }
    },
]

# Get specific frames from a trajectory
# WARNING: This function is time unefficient
def get_trajectory_subset (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    start : int = 0,
    end : int = None,
    step : int = 1
):
    # We need an output trajectory filename
    if not output_trajectory_filename:
        raise SystemExit('Missing output trajectory filename')

    # End must be grater than start
    if end != None and end < start:
        raise SystemExit('End frame must be posterior to start frame')

    # Load the trajectory frame by frame and get only the desired frames
    if input_structure_filename:
        trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    else:
        trajectory = mdt.iterload(input_trajectory_filename, chunk=1)
    # Get the first chunk
    reduced_trajectory = None
    frame_count = 0 # This count works only in case the start frame is out of range, for the logs
    for i, chunk in enumerate(trajectory):
        frame_count = i
        if i == start:
            reduced_trajectory = chunk
            break
    # If we have nothing at this point then it means our start is out of the frames range
    if not reduced_trajectory:
        frame_count += 1
        raise SystemExit('The trajectory has ' + str(frame_count) + ' frames so we can not start at frame ' + str(start))
    # Get further chunks
    relative_end = end - start if end != None else None
    for i, chunk in enumerate(trajectory, 1): # Start the count at 1
        if i == relative_end:
            break
        if i % step == 0:
            reduced_trajectory = mdt.join([reduced_trajectory, chunk], check_topology=False)

    # WARNING: This print here is not just a log. DO NOT REMOVE IT
    # WARNING: It fixes an error (ValueError: unitcell angle < 0) which happens sometimes
    # WARNING: It fixes an error (ValueError: Only rectilinear boxes can be saved to mdcrd files) which happens sometimes
    # DANI: Esto es azaroso: a veces funciona y a veces no
    #print(reduced_trajectory)

    # This is necessary sometimes to avoid the following error:
    #     ValueError: Only rectilinear boxes can be saved to mdcrd files
    reduced_trajectory.unitcell_lengths = [0,0,0]
    reduced_trajectory.unitcell_angles = [90,90,90]

    # Write reduced trajectory to output file
    reduced_trajectory.save(output_trajectory_filename)
get_trajectory_subset.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_trajectory_formats
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': mdtraj_supported_structure_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_structure_formats
        }
    }
]

# Get first frame from a trajectory
def get_first_frame (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_frame_filename : str,
):
    get_trajectory_subset(input_structure_filename, input_trajectory_filename, output_frame_filename, 0)
# Get last last from a trajectory
def get_last_frame (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_frame_filename : str,
):
    frame_count = get_frames_count(input_structure_filename, input_trajectory_filename)
    get_trajectory_subset(input_structure_filename, input_trajectory_filename, output_frame_filename, frame_count-1)
# Set function supported formats
single_frame_getter_format_sets = [
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_frame_filename': mdtraj_supported_trajectory_formats
        }
    },
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_frame_filename': { *mdtraj_supported_structure_formats, *mdtraj_supported_trajectory_formats }
        }
    },
]
get_first_frame.format_sets = single_frame_getter_format_sets
get_last_frame.format_sets = single_frame_getter_format_sets

# Split a trajectory which is actually a merge of independent trajectories back to the original pieces
# Run an RMSD analysis to guess where the pieces are
# The cutoff sets the minimum RMSD change to consider it is a different trajectory
def split_merged_trajectories (
    input_structure_filename : str,
    input_trajectory_filename : str,
    sudden_jump_cutoff : float = 0.2,
    output_trajectory_prefix : str = 'split'):
    # Get the input trajectory format
    input_trajectory_format = get_format(input_trajectory_filename)
    # The cutoff must be a negative number since independent trajectories RMSD sudden jumps will be negative
    cutoff = -abs(sudden_jump_cutoff)
    # Load the trajectory
    trajectory = mdt.load(input_trajectory_filename, top=input_structure_filename)
    # Run a RMSD analysis
    rmsd_data = mdt.rmsd(trajectory, trajectory, 0)
    # Find sudden jumps in RMSD values
    sudden_jumps = []
    previous_value = 0
    for i, value in enumerate(rmsd_data):
        diff = value - previous_value
        # In case this is a new trajectory the jump will be negative
        if diff < cutoff:
            print('New trajectory at frame ' + str(i))
            sudden_jumps.append(i)
        previous_value = value
    # In there was no jumps then stop here
    if len(sudden_jumps) == 0:
        print('Apparently there is a single trajectory')
        return
    # Generate a trajectory subset for each cut
    cut_points = [ 0, *sudden_jumps, len(rmsd_data) ]
    for i in range(len(cut_points) -1):
        start = cut_points[i]
        end = cut_points[i+1]
        trajectory_split = trajectory[start:end]
        split_filename = output_trajectory_prefix + '_' + str(i+1) + '.' + input_trajectory_format
        print('Writting from frame ' + str(start) + ' to frame ' + str(end) + ' to "' + split_filename + '"')
        trajectory_split.save(split_filename)
    
# Sort atoms in a trajectory file by sorting its coordinates
# A new atom indices list is expected in order to do the sorting
def sort_trajectory_atoms (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    new_atom_indices : List[int]
):

    # Assert we have input values
    if not input_structure_filename:
        raise SystemExit('ERROR: Missing input structure filename')
    if not input_trajectory_filename:
        raise SystemExit('ERROR: Missing input trajectory filename')
    if not output_trajectory_filename:
        raise SystemExit('ERROR: Missing output trajectory filename')

    # Load the topology, which is used further
    topology = mdt.load_topology(input_structure_filename)

    # Check the new atom indices to match the number of atoms in the topology
    if len(new_atom_indices) != topology.n_atoms:
        raise ValueError('The number of atom indices (' + str(len(new_atom_indices)) + ') does not match the number of atoms in topology(' + str(topology.n_atoms) + ')' )

    # If the input and the output trajectory filenames are equal then we must rename the input filename to not overwrite it
    if input_trajectory_filename == output_trajectory_filename:
        backup_filename = 'backup.' + input_trajectory_filename
        os.rename(input_trajectory_filename, backup_filename)
        input_trajectory_filename = backup_filename

    # If the output trajectory file already exists at this point then we must stop here
    # The raw trjcat implementation will not join things to the end of it
    if os.path.exists(output_trajectory_filename):
        raise SystemExit('The output file already exists and overwrite is not supported for this functionality')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    frame_filename = '.frame.xtc'

    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()

    for f, frame in enumerate(trajectory):
        # Update the current frame log
        print(ERASE_PREVIOUS_LINE)
        print('Frame ' + str(f))

        # Sort coordinates according to new atom indices
        xyz = frame.xyz[0]
        new_xyz = np.empty(shape=xyz.shape)
        for i, atom_index in enumerate(new_atom_indices):
            new_xyz[i] = xyz[atom_index]

        # Export the current frame
        new_frame = mdt.Trajectory([new_xyz],
            topology=topology,
            time=frame.time, # This is necessary for the new mdtraj trajectory to be functional
            unitcell_lengths=frame.unitcell_lengths, # This is necessary for the new mdtraj trajectory to be functional
            unitcell_angles=frame.unitcell_angles # This is necessary for the new mdtraj trajectory to be functional
        )
        new_frame.save(frame_filename)
        # Join current frame to the output trajectory
        merge_xtc_files(output_trajectory_filename, frame_filename)

    # Remove the residual file
    # WARNING: It may not exist if the trajectory has 1 frame
    if os.path.exists(frame_filename):
        os.remove(frame_filename)