from typing import List, Set, Optional

from .formats import get_format, get_format_set_suitable_function
from .gmx_spells import get_first_frame as gmx_get_first_frame
from .mdt_spells import get_first_frame as mdt_get_first_frame
from .mdt_spells import get_last_frame as mdt_get_last_frame

# Set the functions to perform single frame gettering
first_frame_getter_functions = [ gmx_get_first_frame, mdt_get_first_frame ]
last_frame_getter_functions = [ mdt_get_last_frame ]

# Get the first frame from a trajectory
# Return the single frame filename
# WARNING: The structure filename may be None
def get_first_frame (
    input_structure_filename : str,
    input_trajectory_filename : str,
    accepted_output_formats : Optional[Set[str]] = None
) -> str:
    # Get the input formats
    input_structure_format = get_format(input_structure_filename)
    input_trajectory_format = get_format(input_trajectory_filename)
    format_set = {
        'inputs': {
            'input_structure_filename': None if input_structure_format == None else { input_structure_format },
            'input_trajectory_filename': { input_trajectory_format }
        },
        'outputs': {
            'output_frame_filename': accepted_output_formats
        }
    }
    # Get a suitable function to obtain the unique frame
    suitables = next(get_format_set_suitable_function(
        available_functions=first_frame_getter_functions,
        available_request_format_sets=[format_set],
    ), None)
    if not suitables:
        raise SystemExit('There is no first frame getter function which supports the requested formats')
    suitable_function, formats = suitables
    # The output format will be the first common format between the available formats and the function formats
    output_format = list(formats['outputs']['output_frame_filename'])[0]
    # Set the output filename
    output_single_frame_filename = '.single_frame.' + output_format
    suitable_function(input_structure_filename, input_trajectory_filename, output_single_frame_filename)
    return output_single_frame_filename

# Get the last frame from a trajectory
# Return the single frame filename
# WARNING: The structure filename may be None
def get_last_frame (
    input_structure_filename : str,
    input_trajectory_filename : str,
    accepted_output_formats : Optional[Set[str]] = None
) -> str:
    # Get the input formats
    input_structure_format = get_format(input_structure_filename)
    input_trajectory_format = get_format(input_trajectory_filename)
    format_set = {
        'inputs': {
            'input_structure_filename': None if input_structure_format == None else { input_structure_format },
            'input_trajectory_filename': { input_trajectory_format }
        },
        'outputs': {
            'output_frame_filename': accepted_output_formats
        }
    }
    # Get a suitable function to obtain the unique frame
    suitables = next(get_format_set_suitable_function(
        available_functions=last_frame_getter_functions,
        available_request_format_sets=[format_set],
    ), None)
    if not suitables:
        raise SystemExit('There is no last frame getter function which supports the requested formats')
    suitable_function, formats = suitables
    # The output format will be the last common format between the available formats and the function formats
    output_format = list(formats['outputs']['output_frame_filename'])[0]
    # Set the output filename
    output_single_frame_filename = '.single_frame.' + output_format
    suitable_function(input_structure_filename, input_trajectory_filename, output_single_frame_filename)
    return output_single_frame_filename
