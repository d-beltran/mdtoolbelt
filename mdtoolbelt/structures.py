# Main handler of the toolbelt
import os
from typing import Optional, Tuple, List
Coords = Tuple[float, float, float]

import prody

from .selections import Selection

# An atom
class Atom:
    def __init__ (self,
        name : Optional[str] = None,
        element : Optional[str] = None,
        coords : Optional[Coords] = None,
        residue_index : Optional[int] = None,
        chain_index : Optional[int] = None,
        ):
        self.name = name
        self.element = element
        self.coords = coords
        self.residue_index = residue_index
        self.chain_index = chain_index
        # Set variables to store references to other related instances
        # These variables will be set further by the structure
        self.structure = None
        self.residue = None
        self.chain = None

    def __repr__ (self):
        return '<Atom ' + self.name + '>'


# A residue
class Residue:
    def __init__ (self,
        name : Optional[str] = None,
        number : Optional[int] = None,
        icode : Optional[str] = None,
        atom_indices : List[int] = [],
        chain_index : Optional[int] = None,
        ):
        self.name = name
        self.number = number
        self.icode = icode
        self.atom_indices = atom_indices
        self.chain_index = chain_index
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self.structure = None
        self.atoms = None
        self.chain = None

    def __repr__ (self):
        return '<Residue ' + self.name + str(self.number) + (self.icode if self.icode else '') + '>'


# A chain
class Chain:
    def __init__ (self,
        name : Optional[str] = None,
        atom_indices : List[int] = [],
        residue_indices : List[int] = [],
        ):
        self.name = name
        self.atom_indices = atom_indices
        self.residue_indices = residue_indices
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self.structure = None
        self.atoms = None
        self.residues = None

    def __repr__ (self):
        return '<Chain ' + self.name + '>'


# A structure is a group of atoms organized in chains and residues
class Structure:
    def __init__ (self,
        atoms : List['Atom'] = [],
        residues : List['Residue'] = [],
        chains : List['Chain'] = [],
        ):
        self.atoms = atoms
        self.residues = residues
        self.chains = chains
        # Set self as structure on each input atom, residue and chain
        for instance in atoms + residues + chains:
            instance.structure = self
        # Set references between instances
        for atom in atoms:
            atom.residue = residues[atom.residue_index]
            atom.chain = chains[atom.chain_index]
        for residue in residues:
            residue.atoms = [ atoms[index] for index in residue.atom_indices ]
            residue.chain = chains[residue.chain_index]
        for chain in chains:
            chain.atoms = [ atoms[index] for index in chain.atom_indices ]
            chain.residues = [ residues[index] for index in chain.residue_indices ]
        # Set internal variables
        self._prody_topology = None

    def __repr__ (self):
        return '<Structure (' + str(len(self.atoms)) + ' atoms)>'

    # Set the structure from a ProDy topology
    @classmethod
    def from_prody(cls, prody_topology):
        parsed_atoms = []
        parsed_residues = []
        parsed_chains = []
        prody_atoms = list(prody_topology.iterAtoms())
        prody_residues = list(prody_topology.iterResidues())
        prody_chains = list(prody_topology.iterChains())
        # Parse atoms
        for prody_atom in prody_atoms:
            name = prody_atom.getName()
            element = prody_atom.getElement()
            coords = tuple(prody_atom.getCoords())
            residue_index = prody_atom.getResindex()
            chain_index = prody_atom.getChindex()
            parsed_atom = Atom(name=name, element=element, coords=coords, residue_index=residue_index, chain_index=chain_index)
            parsed_atoms.append(parsed_atom)
        # Parse residues
        for prody_residue in prody_residues:
            name = prody_residue.getResname()
            number = prody_residue.getResnum()
            icode = prody_residue.getIcode()
            atom_indices = list(prody_residue.getIndices())
            chain_index = prody_residue.getChindices()[0]
            parsed_residue = Residue(name=name, number=number, icode=icode, atom_indices=atom_indices, chain_index=chain_index)
            parsed_residues.append(parsed_residue)
        # Parse chains
        for prody_chain in prody_chains:
            name = prody_chain.getChid()
            atom_indices = list(prody_chain.getIndices())
            residue_indices = [ residue.getResindex() for residue in prody_chain.iterResidues() ]
            parsed_chain = Chain(name=name, atom_indices=atom_indices, residue_indices=residue_indices)
            parsed_chains.append(parsed_chain)
        return cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)

    # Set the structure from a pdb file
    # Use ProDy to do so
    @classmethod
    def from_pdb_file(cls, pdb_filename : str):
        if not os.path.exists(pdb_filename):
            raise SystemExit('File "' + pdb_filename + '" not found')
        prody_topology = prody.parsePDB(pdb_filename)
        return cls.from_prody(prody_topology)

    # Generate a pdb file with current structure
    def generate_pdb_file(self, pdb_filename : str):
        with open(pdb_filename, "w") as file:
            file.write('REMARK mdtoolbelt dummy pdb file\n')
            for a, atom in enumerate(self.atoms):
                residue = atom.residue
                index = str(a+1).rjust(5)
                name =  ' ' + atom.name.ljust(3) if len(atom.name) < 4 else atom.name
                residue_name = residue.name.ljust(4)
                chain = atom.chain.name.rjust(1)
                residue_number = str(residue.number).rjust(4)
                icode = residue.icode.rjust(1)
                coords = atom.coords
                x_coord, y_coord, z_coord = [ "{:.3f}".format(coord).rjust(8) for coord in coords ]
                occupancy = '1.00' # Just a placeholder
                temp_factor = '0.00' # Just a placeholder
                element = atom.element
                atom_line = ('ATOM  ' + index + ' ' + name + ' ' + residue_name
                    + chain + residue_number + icode + '   ' + x_coord + y_coord + z_coord
                    + '  ' + occupancy + '  ' + temp_factor + ' ' + element).ljust(80) + '\n'
                file.write(atom_line)

    # Get the structure equivalent prody topology
    def get_prody_topology (self):
        # Return the internal perimeter value if it exists
        if self._prody_topology:
            return self._prody_topology
        # If not, generate the prody topology
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        prody_topology = prody.parsePDB(pdb_filename)
        os.remove(pdb_filename)
        self._prody_topology = prody_topology
        return prody_topology

    # The equivalent prody topology
    prody_topology = property(get_prody_topology, None, None, "The structure equivalent prody topology")

    # Select atoms from the structure thus generating an atom indices list
    # Different tools may be used to make the selection:
    # - prody (default)
    def select (self, selection_string : str, logic : str = 'prody') -> Optional['Selection']:
        if logic == 'prody':
            prody_selection = self.prody_topology.select(selection_string)
            if not prody_selection:
                print('WARNING: Empty selection')
                return None
            return Selection.from_prody(prody_selection)
    
    # Create a new structure from the current using a selection to filter atoms
    def filter (self, selection : 'Selection') -> 'Structure':
        if not selection:
            raise SystemExit('No selection was passed')
        new_atoms = []
        new_residues = []
        new_chains = []
        # Get the selected atoms
        for index in selection.atom_indices:
            # Make a copy of the selected atoms in order to not modify the original ones
            original_atom = self.atoms[index]
            new_atom = Atom(
                name=original_atom.name,
                element=original_atom.element,
                coords=original_atom.coords,
                residue_index=original_atom.residue_index,
                chain_index=original_atom.chain_index
            )
            new_atoms.append(new_atom)
        # Find the selected residues
        selected_residue_indices = list(set([ atom.residue_index for atom in new_atoms ]))
        for index in selected_residue_indices:
            original_residue = self.residues[index]
            new_residue = Residue(
                name=original_residue.name,
                number=original_residue.number,
                icode=original_residue.icode,
                atom_indices=original_residue.atom_indices,
                chain_index=original_residue.chain_index
            )
            new_residues.append(new_residue)
        # Find the selected chains
        selected_chain_indices = list(set([ atom.chain_index for atom in new_atoms ]))
        for index in selected_chain_indices:
            original_chain = self.chains[index]
            new_chain = Chain(
                name=original_chain.name,
                atom_indices=original_chain.atom_indices,
                residue_indices=original_chain.residue_indices
            )
            new_chains.append(new_chain)
        # Generate dictionaries with new indexes as keys and previous indexes as values for atoms, residues and chains
        old_atom_indices = {}
        for i, index in enumerate(selection.atom_indices):
            old_atom_indices[index] = i
        old_residue_indices = {}
        for i, index in enumerate(selected_residue_indices):
            old_residue_indices[index] = i
        old_chain_indices = {}
        for i, index in enumerate(selected_chain_indices):
            old_chain_indices[index] = i
        # Finally, reset indices in all instances
        for atom in new_atoms:
            atom.residue_index = old_residue_indices[atom.residue_index]
            atom.chain_index = old_chain_indices[atom.chain_index]
        for residue in new_residues:
            residue.atom_indices = [ old_atom_indices[index] for index in residue.atom_indices ]
            residue.chain_index = old_chain_indices[residue.chain_index]
        for chain in new_chains:
            chain.atom_indices = [ old_atom_indices[index] for index in chain.atom_indices ]
            chain.residue_indices = [ old_residue_indices[index] for index in chain.residue_indices ]
        return Structure(atoms=new_atoms, residues=new_residues, chains=new_chains)