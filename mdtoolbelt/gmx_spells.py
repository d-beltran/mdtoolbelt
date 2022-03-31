import os
from shutil import copyfile
from subprocess import run, PIPE, Popen

from typing import List

from .formats import get_format

# Get the first frame from a trajectory
def get_first_frame (input_structure_filename : str, input_trajectory_filename : str, output_frame_filename : str):
    # Run Gromacs
    if input_structure_filename:
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_structure_filename,
            "-f",
            input_trajectory_filename,
            "-o",
            output_frame_filename,
            "-dump",
            "0",
            "-quiet"
        ], stdin=p.stdout, stderr=PIPE).stderr.decode()
    else:
        logs = run([
            "gmx",
            "trjconv",
            "-f",
            input_trajectory_filename,
            "-o",
            output_first_frame_filename,
            "-dump",
            "0",
            "-quiet"
        ], stderr=PIPE).stderr.decode()
    # If output has not been generated then warn the user
    if not os.path.exists(output_frame_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

# Set function supported formats
get_first_frame.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr', 'pdb', 'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'pdb', 'gro'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'pdb'}
        },
        'outputs': {
            'output_frame_filename': {'pdb', 'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'gro'}
        },
        'outputs': {
            'output_frame_filename': {'gro', 'xtc', 'trr'}
        }
    }
]

# Get the structure of a tpr file using the first frame getter function
def get_tpr_structure (input_structure_filename : str, input_trajectory_filename : str, output_structure_filename : str):
    get_first_frame(input_structure_filename, input_trajectory_filename, output_structure_filename)
get_tpr_structure.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr', 'pdb', 'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_structure_filename': {'pdb', 'gro'}
        }
    }
]

# Get gromacs supported trajectories merged and converted to a different format
def merge_and_convert_trajectories (input_trajectory_filenames : List[str], output_trajectory_filename : str):
    # Get trajectory formats
    sample_trajectory = input_trajectory_filenames[0]
    input_trajectories_format = get_format(sample_trajectory)
    output_trajectory_format = get_format(output_trajectory_filename)
    auxiliar_single_trajectory_filename = '.single_trajectory.' + input_trajectories_format
    # If we have multiple trajectories then join them
    if len(input_trajectory_filenames) > 1:
        single_trajectory_filename = auxiliar_single_trajectory_filename
        logs = run([
            "gmx",
            "trjcat",
            "-f",
            *input_trajectory_filenames,
            "-o",
            single_trajectory_filename,
            "-quiet"
        ], stderr=PIPE).stderr.decode()
        # If output has not been generated then warn the user
        if not os.path.exists(single_trajectory_filename):
            print(logs)
            raise SystemExit('Something went wrong with Gromacs')
    else:
        single_trajectory_filename = sample_trajectory
    # In case input and output formats are different we must convert the trajectory
    if input_trajectories_format != output_trajectory_format:
        logs = run([
            "gmx",
            "trjconv",
            "-f",
            single_trajectory_filename,
            "-o",
            output_trajectory_filename,
            "-quiet"
        ], stderr=PIPE).stderr.decode()
        # If output has not been generated then warn the user
        if not os.path.exists(output_trajectory_filename):
            print(logs)
            raise SystemExit('Something went wrong with Gromacs')
    else:
        copyfile(single_trajectory_filename, output_trajectory_filename)
    # Remove residual files
    if os.path.exists(auxiliar_single_trajectory_filename):
        os.remove(auxiliar_single_trajectory_filename)

merge_and_convert_trajectories.format_sets = [
    {
        'inputs': {
            'input_trajectory_filenames': {'xtc', 'trr'}
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_trajectory_filenames': {'pdb'}
        },
        'outputs': {
            'output_trajectory_filename': {'pdb', 'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_trajectory_filenames': {'gro'}
        },
        'outputs': {
            'output_trajectory_filename': {'gro', 'xtc', 'trr'}
        }
    }
]