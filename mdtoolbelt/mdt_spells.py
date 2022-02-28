import os
from subprocess import run, PIPE
import mdtraj as md

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
    trajectory = md.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
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
            reduced_trajectory = md.join([reduced_trajectory, chunk], check_topology=False)

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