from os.path import exists
from .formats import get_format
from .gmx_spells import filter_structure, filter_trajectory
from .structures import Structure

# Note tha this function only accepts pdb/gro structures and xtc/trr trajectories
accepted_structure_formats = ['pdb']
accepted_trajectory_formats = ['xtc', 'trr']

# Filter both structure and trajectory
# Note that this function is not format-smart
def filter_atoms (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_structure_filename : str,
    output_trajectory_filename : str,
    selection_string : str,
    selection_syntax : str = 'vmd'
):  

    # Check formats are as expected
    # Check also input files exist
    if not input_structure_filename:
        raise SystemExit('Missing input structure filename')
    input_structure_format = get_format(input_structure_filename)
    if input_structure_format not in accepted_structure_formats:
        raise SystemExit('Not valid input structure format (' + input_structure_format + '). Accepted structure formats: ' + ','.join(accepted_structure_formats))
    if not exists(input_structure_filename):
        raise SystemExit('Missing input file ' + input_structure_filename)
    if input_trajectory_filename:
        input_trajectory_format = get_format(input_trajectory_filename)
        if input_trajectory_format not in accepted_trajectory_formats:
            raise SystemExit('Not valid input trajectory format (' + input_trajectory_format + '). Accepted trajectory formats: ' + ','.join(accepted_trajectory_formats))
        if not exists(input_trajectory_filename):
            raise SystemExit('Missing input file ' + input_trajectory_filename)
    if output_structure_filename:
        output_structure_format = get_format(output_structure_filename)
        if output_structure_format not in accepted_structure_formats:
            raise SystemExit('Not valid output structure format (' + output_structure_format + '). Accepted structure formats: ' + ','.join(accepted_structure_formats))
    if output_trajectory_filename:
        output_trajectory_format = get_format(output_trajectory_filename)
        if output_trajectory_format not in accepted_trajectory_formats:
            raise SystemExit('Not valid output trajectory format (' + output_trajectory_format + '). Accepted trajectory formats: ' + ','.join(accepted_trajectory_formats))
    # Check we are not missing additional inputs
    if not selection_string:
        raise SystemExit('Missing input selection string')
    if not selection_syntax:
        raise SystemExit('Missing input selection syntax')

    # Parse the selection
    structure = Structure.from_pdb_file(input_structure_filename)
    selection = structure.select(selection_string, syntax=selection_syntax)

    # If the selection is empty then war the user
    if not selection:
        raise SystemExit('Selection ' + selection_string + ' is empty')

    # Run the filters
    if output_structure_filename:
        filter_structure(input_structure_filename, output_structure_filename, selection)
    if input_trajectory_filename and output_trajectory_filename:
        filter_trajectory(input_structure_filename, input_trajectory_filename, output_trajectory_filename, selection)