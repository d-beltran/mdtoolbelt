from .formats import get_format, get_format_set_suitable_function
from .pyt_spells import get_frames_count as pyt_get_frames_count

# Set the functions to perform trajectory frame count
frame_count_functions = [ pyt_get_frames_count ]

# Get the trajectory frame count
# WARNING: The structure filename may be None
def get_frames_count (
    input_structure_filename : str,
    input_trajectory_filename : str,
) -> int:
    # Get the input formats
    input_structure_format = get_format(input_structure_filename)
    input_trajectory_format = get_format(input_trajectory_filename)
    format_set = {
        'inputs': {
            'input_structure_filename': None if input_structure_format == None else { input_structure_format },
            'input_trajectory_filename': { input_trajectory_format }
        }
    }
    # Get a suitable function
    suitables = next(get_format_set_suitable_function(
        available_functions=frame_count_functions,
        available_request_format_sets=[format_set],
    ), None)
    if not suitables:
        raise SystemExit('There is no frame count function which supports the requested formats')
    suitable_function, formats = suitables
    # Run the suitable function
    frame_count = suitable_function(input_structure_filename, input_trajectory_filename)
    return frame_count