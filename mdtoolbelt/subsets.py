from inspect import getfullargspec

from .formats import get_format, get_format_set_suitable_function
from .gmx_spells import get_trajectory_subset as gmx_get_trajectory_subset
from .mdt_spells import get_trajectory_subset as mdt_get_trajectory_subset

# Set the functions to perform single frame gettering
subset_functions = [ gmx_get_trajectory_subset, mdt_get_trajectory_subset ]

# Get specific frames from a trajectory
def get_trajectory_subset (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    start : int = 0,
    end : int = None,
    step : int = 1
):
    # Get the input formats
    input_structure_format = get_format(input_structure_filename)
    input_trajectory_format = get_format(input_trajectory_filename)
    output_trajectory_format = get_format(output_trajectory_filename)
    format_set = {
        'inputs': {
            'input_structure_filename': None if input_structure_format == None else { input_structure_format },
            'input_trajectory_filename': { input_trajectory_format }
        },
        'outputs': {
            'output_trajectory_filename': { output_trajectory_format }
        }
    }
    # Get a suitable function to do the job
    suitables = next(get_format_set_suitable_function(
        available_functions=subset_functions,
        available_request_format_sets=[format_set],
    ), None)
    if not suitables:
        raise SystemExit('There is no subset function which supports the requested formats')
    suitable_function, formats = suitables
    # Call the subset function
    # Check if the subset function requires the input structure argument before
    suitable_function_keywords = getfullargspec(suitable_function)[0]
    required_structure = 'input_structure_filename' in suitable_function_keywords
    if required_structure:
        suitable_function(input_structure_filename, input_trajectory_filename, output_trajectory_filename, start, end, step)
    else:
        suitable_function(input_trajectory_filename, output_trajectory_filename, start, end, step)