from typing import Optional, List, Tuple, Callable

# Get a filename format
def get_format (filename : str) -> str:
    return filename.split('.')[-1]

# Find a function which is suitable for any of the available request "format sets"**
# All functions are checked for each request format set before jumping to another and they are evaluated in order
# A function and new generated format set with formats in common are returned
# None is returned when there is no suitable function
# WARNING: Available functions must have the 'format_sets' property
# ** Format sets are dictionaries which specify input and output formats
# Consider function format sets as 'required input' and 'available output'
# Consider request format sets as 'available input' and 'required output'
# Both inputs and ouputs are dictionaries where keys are function keywords and values are sets of supported formats
# Alternatively, a keyword may have None as value to represent unnecessary requirements or missing availabilities
# An example is shown below:
# {
#     'inputs': {
#         'input_structure_filename': {'tpr'},
#         'input_trajectory_filename': {'xtc', 'trr'},
#     },
#     'outputs': {
#         'output_trajectory_filename': {'pdb', 'gro'}
#     },
# }
# Functions may have multiple format sets since different input formats may lead to different output formats
def get_format_set_suitable_function (
    available_functions : List[Callable],
    available_request_format_sets : List[dict],
) -> Optional[ Tuple[ Callable, dict ] ]:
    # Try with each request format set
    for request_format_set in available_request_format_sets:
        # Search functions to match formats for every keyword
        for function in available_functions:
            # Test every function format set independently
            for function_format_set in function.format_sets:
                # Get common values for every keyword format set
                common_format_set = {'inputs':{},'outputs':{}}
                missing_common_format = False
                # Check format keys are compatible
                if not check_format_sets_compability(request_format_set, function_format_set):
                    raise ValueError('Format keys are not compatible with function ' + str(function.__name__))
                # Check the function inputs to be fulfilled by the request inputs
                required_input_keywords = function_format_set['inputs'].keys()
                for keyword in required_input_keywords:
                    required_formats = function_format_set['inputs'][keyword]
                    # If this format set does not need an input we set None for this keyword
                    if required_formats == None:
                        common_format_set[keyword] = None
                        continue
                    available_formats = request_format_set['inputs'][keyword]
                    # If this format set is missing an input we are done
                    if available_formats == None:
                        missing_common_format = True
                        break
                    common_formats = available_formats.intersection(required_formats)
                    # If the common formats set is empty we are done
                    if not common_formats:
                        missing_common_format = True
                        break
                    common_format_set['inputs'][keyword] = common_formats
                # If any of the common format sets was empty it means formats do not match
                if missing_common_format:
                    continue
                # Check the request outputs to be fulfilled by the function outputs
                required_outpt_keywords = request_format_set['outputs'].keys()
                for keyword in required_outpt_keywords:
                    required_formats = request_format_set['outputs'][keyword]
                    # If this format set does not need an input we set None for this keyword
                    if required_formats == None:
                        common_format_set[keyword] = None
                        continue
                    available_formats = function_format_set['outputs'][keyword]
                    # If this format set is missing an input we are done
                    if available_formats == None:
                        missing_common_format = True
                        break
                    common_formats = available_formats.intersection(required_formats)
                    # If the common formats set is empty we are done
                    if not common_formats:
                        missing_common_format = True
                        break
                    common_format_set['outputs'][keyword] = common_formats
                # If any of the common format sets was empty it means formats do not match
                if missing_common_format:
                    continue
                # Otherwise we have the function
                return function, common_format_set
                

# Check two format sets to be compatible
# Both function and request format sets must match in their requirements
# i.e. all function format set input keywords must be included in request format set input keywords
# i.e. all request format set output keywords must be included in function format set output keywords
def check_format_sets_compability (request_format_set : dict, function_format_set : dict) -> bool:
    # Check the function inputs keyowrds to exist in the request input keywords
    required_input_keywords = function_format_set['inputs'].keys()
    available_input_keywords = request_format_set['inputs'].keys()
    for keyword in required_input_keywords:
        if keyword not in available_input_keywords:
            print('ERROR: Missing ' + keyword + ' keyword')
            return False
    # Check the request output keyowrds to exist in the function output keywords
    required_outpt_keywords = request_format_set['outputs'].keys()
    available_outpt_keywords = function_format_set['outputs'].keys()
    for keyword in required_outpt_keywords:
        if keyword not in available_outpt_keywords:
            print('ERROR: Missing ' + keyword + ' keyword')
            return False
    return True




# Structure file formats

def is_pdb (filename : str) -> bool:
    return filename[-4:] == '.pdb'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

def is_tpr (filename : str) -> bool:
    return filename[-4:] == '.tpr'

def is_gro (filename : str) -> bool:
    return filename[-4:] == '.gro'

def is_prmtop (filename : str) -> bool:
    return filename[-7:] == '.prmtop'

def is_top (filename : str) -> bool:
    return filename[-4:] == '.top'

# Trajectory file formats

def is_xtc (filename : str) -> bool:
    return filename[-4:] == '.xtc'

def is_dcd (filename : str) -> bool:
    return filename[-4:] == '.dcd'

def is_netcdf (filename : str) -> bool:
    return filename[-3:] == '.nc'

def are_xtc (filenames : list) -> bool:
    return all([ is_xtc(filename) for filename in filenames ])

def are_dcd (filenames : list) -> bool:
    return all([ is_dcd(filename) for filename in filenames ])

def are_netcdf (filenames : list) -> bool:
    return all([ is_netcdf(filename) for filename in filenames ])

# Extra formats logic

# Check if a file may be read by pytraj according to its format
def is_pytraj_supported (filename : str) -> bool:
    return is_prmtop(filename) or is_top(filename) or is_psf(filename)

# From GitHub:
# ParmFormatDict = {
#     "AMBERPARM": AMBERPARM,
#     "PDBFILE": PDBFILEPARM,
#     "MOL2FILE": MOL2FILEPARM,
#     "CHARMMPSF": CHARMMPSF,
#     "CIFFILE": CIFFILE,
#     "GMXTOP": GMXTOP,
#     "SDFFILE": SDFFILE,
#     "TINKER": TINKERPARM,
#     "UNKNOWN_PARM": UNKNOWN_PARM,
# }

# Get the pytraj format key for the write_parm function for a specific file according to its format
def get_pytraj_parm_format (filename : str) -> str:
    if is_prmtop(filename):
        return 'AMBERPARM'
    if is_psf(filename):
        return 'CHARMMPSF'
    if is_top(filename):
        return 'GMXTOP'
    if is_pdb(filename):
        return 'PDBFILE'
    raise ValueError('The file ' + filename + ' format is not supported')