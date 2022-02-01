from typing import List

# A selection is a list of atom indices from a structure
class Selection:

    def __init__ (self, atom_indices : List[int]):
        self.atom_indices = atom_indices

    def __repr__ (self):
        return '<Selection (' + str(len(self.atom_indices)) + ' atoms)>'

    @classmethod
    def from_prody (cls, prody_selection):
        indexes = [ atom.getIndex() for atom in prody_selection.iterAtoms() ]
        return cls(indexes)

    # Generate a vmd selection in tcl format
    def to_vmd (self) -> str:
        return 'index ' + ' '.join([ str(index) for index in self.atom_indices ])

    # Generate the file content of gromacs ndx file
    def to_ndx (self, selection_name : str = 'Selection') -> str:
        # Add a header 
        content = '[ ' + selection_name + ' ]\n'
        count = 0
        for index in self.atom_indices:
            # Add a breakline each 15 indices
            count += 1
            if count == 15:
                content += '\n'
                count = 0
            # Add a space between indices
            # Atom indices go from 0 to n-1
            # Add +1 to the index since gromacs counts from 1 to n
            content += str(index + 1) + ' '
        return content