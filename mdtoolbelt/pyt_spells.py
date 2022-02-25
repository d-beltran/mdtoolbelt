import pytraj as pyt

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