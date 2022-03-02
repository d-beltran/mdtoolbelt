import os
from subprocess import run, PIPE
import mdtraj as mdt

from .formats import get_format

# Multiple files may be selected with bash syntax (e.g. *.dcd)
# Tested supported input formats are .dcd
# Tested supported output formats are .xtc
def merge_and_convert_traj (
    input_filenames : list,
    output_filename : str
    ):

    # Run MDtraj
    logs = run([
        "mdconvert",
        "-o",
        output_filename,
        *input_filenames,
    ], stderr=PIPE).stderr.decode()
    # If output has not been generated then warn the user
    if not os.path.exists(output_filename):
        print(logs)
        raise SystemExit('Something went wrong with MDTraj')

def convert_traj (input_trajectory_filename : str, output_trajectory_filename : str):
    merge_and_convert_traj([input_trajectory_filename], output_trajectory_filename)
convert_traj.format_sets = [
    {
        'inputs': {
            'input_trajectory_filename': {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos'}
        },
        'outputs': {
            'output_trajectory_filename': {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos'}
        }
    },
]

# Get specific frames from a trajectory
def get_trajectory_subset (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    start : int = 0,
    end : int = 0,
    step : int = 1
):
    # In case no end is passed return the only the start frame
    if not end:
        end = start + 1

    # Load the trajectory frame by frame and get only the desired frames
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    # Get the first chunk
    reduced_trajectory = None
    for i, chunk in enumerate(trajectory):
        if i == start:
            reduced_trajectory = chunk
            break
    # Get further chunks
    relative_end = end - start
    for i, chunk in enumerate(trajectory, 1): # Start the count at 1
        if i == relative_end:
            break
        if i % step == 0:
            reduced_trajectory = mdt.join([reduced_trajectory, chunk], check_topology=False)

    # Write reduced trajectory to output file
    reduced_trajectory.save(output_trajectory_filename)
get_trajectory_subset.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'pdb', 'h5'},
            'input_trajectory_filename': {'dcd', 'xtc', 'trr', 'nc', 'binpos'}
        },
        'outputs': {
            'output_trajectory_filename': {'dcd', 'xtc', 'trr', 'nc', 'binpos'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'pdb', 'h5'}
        },
        'outputs': {
            'output_trajectory_filename': {'pdb', 'h5'}
        }
    }
]

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
    