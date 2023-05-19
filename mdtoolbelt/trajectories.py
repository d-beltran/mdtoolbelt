from .pyt_spells import get_pytraj_trajectory

# DANI: Esto no se usa de momento

# A trajectory is a group of coordinates
# This class is based in a xtc trajectory file stored in disk and read in a framed way
class Trajectory:
    def __init__ (self, trajectory_filename : str, structure_filename : str = None):
        self.trajectory_filename = trajectory_filename
        self.structure_filename = structure_filename
        self._frame_count = None
        self._pytraj_trajectory = None

    def __repr__ (self):
        return '<Trajectory (' + str(len(self.atoms)) + ' atoms)>'

    # The number of frames in the trajectory (read only)
    def get_frame_count (self) -> int:
        # Return the stored value, if exists
        if self._frame_count != None:
            return self._frame_count
        # If not, we must count the number of frames in the trajectory and store it
        # Note that counting frames may be a long process for huge trajectories
        self._frame_count = self.pytraj_trajectory.n_frames
        return self._frame_count
    frame_count = property(get_frame_count, None, None, "The number of frames in the trajectory (read only)")

    # Trajectory in pytraj format (read only)
    def get_pytraj_trajectory (self):
        # Return the stored value, if exists
        if self._pytraj_trajectory:
            return self._pytraj_trajectory
        # If not, set the pytraj iterloader
        self._pytraj_trajectory = get_pytraj_trajectory(structure_filename, trajectory_filename)
        return self._pytraj_trajectory
    pytraj_trajectory = property(get_pytraj_trajectory, None, None, "Trajectory in pytraj format (read only)")