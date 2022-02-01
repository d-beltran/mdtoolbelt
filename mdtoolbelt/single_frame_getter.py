from typing import List, Optional

from .formats import get_format, get_format_set_suitable_function
from .gmx_spells import gmx_get_first_frame

# Set the functions to performe first frame gettering

# Get the first frame from a trajectory
# Return the single frame filename
first_frame_getter_functions = [ gmx_get_first_frame ]
def get_first_frame ( input_trajectory_filename : str, accepted_output_formats : Optional[List[str]] = None ) -> str:
    # Get the input trajecotry format
    input_format = get_format(input_trajectory_filename)
    format_set = {
        'input_structure_filename': None,
        'input_trajectory_filename': [input_format],
        'output_frame_filename': accepted_output_formats
    }
    # Get a suitable function to obtain the unique frame
    suitables = get_format_set_suitable_function(
        available_functions=first_frame_getter_functions,
        available_request_format_sets=[format_set],
    )
    if not suitables:
        raise SystemExit('There is no first frame getter function which supports the requested formats')
    suitable_function, formats = suitables
    # The output format will be the first common format between the available formats and the function formats
    output_format = formats['outputs']['output_frame_filename'][0]
    # Set the output filename
    output_single_frame_filename = '.single_frame.' + output_format
    suitable_function(input_trajectory_filename, output_single_frame_filename)
    return output_single_frame_filename
