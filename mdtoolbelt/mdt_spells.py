import os
from subprocess import run, PIPE

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