import pytraj as pyt

# Set pytraj supported formats
pytraj_supported_structure_formats = {'prmtop', 'pdb', 'parm7', 'mol2', 'psf', 'cif', 'top', 'sdf'}
pytraj_supported_trajectory_formats = {'xtc', 'trr', 'crd', 'nc', 'dcd'}

# Get the trajectory frames number using pytraj
def get_frames_count (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:
    
    # Load the trajectory from pytraj
    pyt_trajectory = pyt.iterload(
        input_trajectory_filename,
        input_topology_filename)

    # Return the frames number
    return pyt_trajectory.n_frames
# Set function supported formats
get_frames_count.format_sets = [
    {
        'inputs': {
            'input_structure_filename': pytraj_supported_structure_formats,
            'input_trajectory_filename': pytraj_supported_trajectory_formats
        }
    }
]