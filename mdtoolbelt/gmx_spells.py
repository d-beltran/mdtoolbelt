import os
from subprocess import run, PIPE, Popen

from .formats import get_format

# Get the first frame from a trajectory
def gmx_get_first_frame (input_structure_filename : str, input_trajectory_filename : str, output_frame_filename : str):
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
gmx_get_first_frame.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'pdb', 'gro'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': {'pdb'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'pdb'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': {'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'gro'}
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
    gmx_get_first_frame(input_structure_filename, input_trajectory_filename, output_structure_filename)
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