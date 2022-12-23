# Main handler of the toolbelt
import os
import math
from bisect import bisect
from typing import Optional, Union, Tuple, List, Generator
Coords = Tuple[float, float, float]

import prody
import pytraj

from .formats import get_format
from .selections import Selection
from .vmd_spells import get_vmd_selection_atom_indices, get_covalent_bonds
from .mdt_spells import sort_trajectory_atoms
from .utils import residue_name_to_letter

# DANI: git se ha rallado

# Set all available chains according to pdb standards
available_caps = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
    'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

# An atom
class Atom:
    def __init__ (self,
        name : Optional[str] = None,
        element : Optional[str] = None,
        coords : Optional[Coords] = None,
        ):
        self.name = name
        self.element = element
        self.coords = coords
        # Set variables to store references to other related instances
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._residue_index = None

    def __repr__ (self):
        return '<Atom ' + self.name + '>'

    def __eq__ (self, other):
        return self._residue_index == other._residue_index and self.name == other.name

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    index = property(get_index, None, None, "The residue index according to parent structure residues (read only)")

    # The atom residue index according to parent structure residues
    # If residue index is set then make changes in all the structure to make this change coherent
    def get_residue_index (self) -> int:
        return self._residue_index
    def set_residue_index (self, new_residue_index : int):
        # If there is not strucutre yet it means the residue is beeing set before the structure
        # We just save the residue index and wait for the structure to be set
        if not self.structure:
            self._residue_index = new_residue_index
            return
        # Relational indices are updated through a top-down hierarchy
        # Affected residues are the ones to update this atom internal residue index
        current_residue = self.residue
        current_residue.remove_atom(self)
        new_residue = self.structure.residues[new_residue_index]
        new_residue.add_atom(self)
    residue_index = property(get_residue_index, set_residue_index, None, "The atom residue index according to parent structure residues")

    # The atom residue
    # If residue is set then make changes in all the structure to make this change coherent
    def get_residue (self) -> Optional['Residue']:
        # If there is not strucutre yet it means the atom is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        # It also may happend that the residue is missing because the atom has been set rawly
        if not self.structure or self.residue_index == None:
            return None
        # Get the residue in the structure according to the residue index
        return self.structure.residues[self.residue_index]
    def set_residue (self, new_residue : 'Residue'):
        # Find the new residue index and set it as the atom residue index
        # Note that the residue must be set in the structure already
        new_residue_index = new_residue.index
        if new_residue_index == None:
            raise ValueError('Residue ' + str(new_residue) + ' is not set in the structure')
        self.set_residue_index(new_residue_index)
    residue = property(get_residue, set_residue, None, "The atom residue")

    # The atom chain index according to parent structure chains (read only)
    # In order to change the chain index it must be changed in the atom residue
    def get_chain_index (self) -> Optional[int]:
        # The residue may be missing if the atom has been set rawly
        if not self.residue:
            return None
        return self.residue.chain_index
    chain_index = property(get_chain_index, None, None, "The atom chain index according to parent structure chains (read only)")

    # The atom chain (read only)
    # In order to change the chain it must be changed in the atom residue
    def get_chain (self) -> Optional['Chain']:
        # If there is not strucutre yet it means the atom is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        # It also may happend that the chain is missing because the atom has been set rawly
        if not self.structure or self.chain_index == None:
            return None
        # Get the chain in the structure according to the chain index
        return self.structure.chains[self.chain_index]
    chain = property(get_chain, None, None, "The atom chain (read only)")

    # Generate a selection for this atom
    def get_selection (self) -> 'Selection':
        return Selection([self.index])

    # Get indices of other atoms in the structure which are covalently bonded to this atom
    def get_bonds (self) -> Optional[ List[int] ]:
        if not self.structure:
            raise ValueError('The atom has not a structure defined')
        if self.index == None:
            raise ValueError('The atom has not an index defined')
        return self.structure.bonds[self.index]

    # Make a copy of the current atom
    def copy (self) -> 'Atom':
        atom_copy = Atom(self.name, self.element, self.coords)
        atom_copy._structure = self._structure
        atom_copy._index = self._index
        atom_copy._residue_index = self._residue_index
        return atom_copy

# A residue
class Residue:
    def __init__ (self,
        name : Optional[str] = None,
        number : Optional[int] = None,
        icode : Optional[str] = None,
        ):
        self.name = name
        self.number = number
        self.icode = icode
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._atom_indices = []
        self._chain_index = None

    def __repr__ (self):
        return '<Residue ' + self.name + str(self.number) + (self.icode if self.icode else '') + '>'

    def __eq__ (self, other):
        if type(self) != type(other):
            return False
        return (
            self._chain_index == other._chain_index and
            #self.name == other.name and
            self.number == other.number and
            self.icode == other.icode
        )

    def __hash__ (self):
        # WARNING: This is susceptible to duplicated residues
        return hash((self._chain_index, self.number, self.icode))
        # WARNING: This is not susceptible to duplicated residues
        #return hash(tuple(self._atom_indices))

    def same_inputs_as (self, other) -> bool:
        return (
            self.name == other.name and
            self.number == other.number and
            self.icode == other.icode
        )

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    def set_index (self, index):
        # Update residue atoms
        for atom in self.atoms:
            atom._residue_index = index
        # Update residue chain
        current_index = self._index
        chain = self.chain
        residue_chain_index = chain._residue_indices.index(current_index)
        chain._residue_indices[residue_chain_index] = index
        # Finally update self index
        self._index = index
    index = property(get_index, set_index, None, "The residue index according to parent structure residues (read only)")

    # The atom indices according to parent structure atoms for atoms in this residue
    # If atom indices are set then make changes in all the structure to make this change coherent
    def get_atom_indices (self) -> List[int]:
        return self._atom_indices
    def set_atom_indices (self, new_atom_indices : List[int]):
        # If there is not strucutre yet it means the residue is beeing set before the structure
        # We just save atom indices and wait for the structure to be set
        if not self.structure:
            self._atom_indices = new_atom_indices
            return
        # Update the current atoms
        for atom in self.atoms:
            atom._residue_index = None
        # Update the new atoms
        for index in new_atom_indices:
            atom = self.structure.atoms[index]
            atom._residue_index = self.index
        # Now new indices are coherent and thus we can save them
        self._atom_indices = new_atom_indices
    atom_indices = property(get_atom_indices, set_atom_indices, None, "The atom indices according to parent structure atoms for atoms in this residue")

    # The atoms in this residue
    # If atoms are set then make changes in all the structure to make this change coherent
    def get_atoms (self) -> List['Atom']:
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # In this case it is not possible to get related atoms in the structure
        if not self.structure:
            return []
        # Get atoms in the structure according to atom indices
        atoms = self.structure.atoms
        return [ atoms[atom_index] for atom_index in self.atom_indices ]
    def set_atoms (self, new_atoms : List['Atom']):
        # Find indices for new atoms and set their indices as the new atom indices
        # Note that atoms must be set in the structure already
        new_atom_indices = []
        for new_atom in new_atoms:
            new_atom_index = new_atom.index
            if new_atom_index == None:
                raise ValueError('Atom ' + str(new_atom) + ' is not set in the structure')
            new_atom_indices.append(new_atom_index)
        self.set_atom_indices(new_atom_indices)
    atoms = property(get_atoms, set_atoms, None, "The atoms in this residue")

    # Add an atom to the residue
    def add_atom (self, new_atom : 'Atom'):
        # Insert the new atom index in the list of atom indices keeping the order
        new_atom_index = new_atom.index
        sorted_atom_index = bisect(self.atom_indices, new_atom_index)
        self.atom_indices.insert(sorted_atom_index, new_atom_index)
        # Update the atom internal index
        new_atom._residue_index = self.index

    # Remove an atom from the residue
    def remove_atom (self, current_atom : 'Atom'):
        # Remove the current atom index from the atom indices list
        self.atom_indices.remove(current_atom.index) # This index MUST be in the list
        # Update the atom internal index
        current_atom._residue_index = None

    # The residue chain index according to parent structure chains
    # If chain index is set then make changes in all the structure to make this change coherent
    def get_chain_index (self) -> int:
        return self._chain_index
    def set_chain_index (self, new_chain_index : int):
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # We just save the chain index and wait for the structure to be set
        if not self.structure:
            self._chain_index = new_chain_index
            return
        # Relational indices are updated through a top-down hierarchy
        # Affected chains are the ones to update this residue internal chain index
        # WARNING: It is critical to find both current and new chains before removing/adding residues
        # WARNING: It may happend that we remove the last residue in the current chain and the current chain is purged
        # WARNING: Thus the 'new_chain_index' would be obsolete since the structure.chains list would have changed
        current_chain = self.chain
        new_chain = self.structure.chains[new_chain_index]
        current_chain.remove_residue(self)
        new_chain.add_residue(self)
    chain_index = property(get_chain_index, set_chain_index, None, "The residue chain index according to parent structure chains")

    # The residue chain
    # If chain is set then make changes in all the structure to make this change coherent
    def get_chain (self) -> 'Chain':
        # If there is not strucutre yet it means the residue is beeing set before the structure
        # In this case it is not possible to get related chain in the structure
        if not self.structure:
            return []
        # Get the chain in the structure according to the chain index
        return self.structure.chains[self.chain_index]
    def set_chain (self, new_chain : Union['Chain', str]):
        # In case the chain is just a string we must find/create the corresponding chain
        if type(new_chain) == str:
            letter = new_chain
            # Get the residue structure
            structure = self.structure
            if not structure:
                raise ValueError('Cannot find the corresponding ' + new_chain + ' chain without the structure')
            # Find if the letter belongs to an already existing chain
            new_chain = structure.get_chain_by_name(letter)
            # If the chain does not exist yet then create it
            if not new_chain:
                new_chain = Chain(name=letter)
                structure.set_new_chain(new_chain)
        # Find the new chain index and set it as the residue chain index
        # Note that the chain must be set in the structure already
        new_chain_index = new_chain.index
        if new_chain_index == None:
            raise ValueError('Chain ' + str(new_chain) + ' is not set in the structure')
        self.set_chain_index(new_chain_index)
    chain = property(get_chain, set_chain, None, "The residue chain")

    # Generate a selection for this residue
    def get_selection (self) -> 'Selection':
        return Selection(self.atom_indices)

    # Make a copy of the current residue
    def copy (self) -> 'Residue':
        residue_copy = Residue(self.name, self.number, self.icode)
        residue_copy._structure = self._structure
        residue_copy._index = self._index
        residue_copy._atom_indices = self._atom_indices
        residue_copy._chain_index = self._chain_index
        return residue_copy

# A chain
class Chain:
    def __init__ (self,
        name : Optional[str] = None,
        ):
        self.name = name
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._residue_indices = []

    def __repr__ (self):
        return '<Chain ' + self.name + '>'

    def __eq__ (self, other):
        return self.name == other.name

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    # When the index is set all residues are updated with the next chain index
    def set_index (self, index : int):
        for residue in self.residues:
            residue._chain_index = index
        self._index = index
    index = property(get_index, set_index, None, "The residue index according to parent structure residues (read only)")

    # The residue indices according to parent structure residues for residues in this chain
    # If residue indices are set then make changes in all the structure to make this change coherent
    def get_residue_indices (self) -> List[int]:
        return self._residue_indices
    def set_residue_indices (self, new_residue_indices : List[int]):
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # We just save residue indices and wait for the structure to be set
        if not self.structure:
            self._residue_indices = new_residue_indices
            return
        # Update the current residues
        for residue in self.residues:
            residue._chain_index = None
        # Update the new residues
        for index in new_residue_indices:
            residue = self.structure.residues[index]
            residue._chain_index = self.index
        # In case the new residue indices list is empty this chain must be removed from its structure
        if len(new_residue_indices) == 0:
            self.structure.purge_chain(self)
        # Now new indices are coherent and thus we can save them
        self._residue_indices = new_residue_indices
    residue_indices = property(get_residue_indices, set_residue_indices, None, "The residue indices according to parent structure residues for residues in this residue")

    # The residues in this chain
    # If residues are set then make changes in all the structure to make this change coherent
    def get_residues (self) -> List['Residue']:
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        if not self.structure:
            return []
        # Get residues in the structure according to residue indices
        residues = self.structure.residues
        return [ residues[residue_index] for residue_index in self.residue_indices ]
    def set_residues (self, new_residues : List['Residue']):
        # Find indices for new residues and set their indices as the new residue indices
        # Note that residues must be set in the structure already
        new_residue_indices = []
        for new_residue in new_residues:
            new_residue_index = new_residue.index
            if new_residue_index == None:
                raise ValueError('Residue ' + str(new_residue) + ' is not set in the structure')
            new_residue_indices.append(new_residue_index)
        self.set_residue_indices(new_residue_indices)
    residues = property(get_residues, set_residues, None, "The residues in this chain")

    # Add a residue to the chain
    def add_residue (self, residue : 'Residue'):
        # Insert the new residue index in the list of residue indices keeping the order
        sorted_residue_index = bisect(self.residue_indices, residue.index)
        self.residue_indices.insert(sorted_residue_index, residue.index)
        # Update the residue internal chain index
        residue._chain_index = self.index

    # Remove a residue from the chain
    # WARNING: Note that this function does not trigger the set_residue_indices
    def remove_residue (self, residue : 'Residue'):
        self.residue_indices.remove(residue.index) # This index MUST be in the list
        # If we removed the last index then this chain must be removed from its structure
        if len(self.residue_indices) == 0 and self.structure:
            self.structure.purge_chain(self)
        # Update the residue internal chain index
        residue._chain_index = None

    # Atom indices for all atoms in the chain (read only)
    # In order to change atom indices they must be changed in their corresponding residues
    def get_atom_indices (self) -> List[int]:
        return sum([ residue.atom_indices for residue in self.residues ], [])
    atom_indices = property(get_atom_indices, None, None, "Atom indices for all atoms in the chain (read only)")

    # Atoms in the chain (read only)
    # In order to change atoms they must be changed in their corresponding residues
    def get_atoms (self) -> List[int]:
        return sum([ residue.atoms for residue in self.residues ], [])
    atoms = property(get_atoms, None, None, "Atoms in the chain (read only)")

    # Get the residues sequence in one-letter code
    def get_sequence (self) -> str:
        return ''.join([ residue_name_to_letter(residue.name) for residue in self.residues ])

    # Generate a selection for this chain
    def get_selection (self) -> 'Selection':
        return Selection(self.atom_indices)

    # Make a copy of the current chain
    def copy (self) -> 'Chain':
        chain_copy = Chain(self.name)
        chain_copy._structure = self._structure
        chain_copy._index = self._index
        chain_copy.residue_indices = self.residue_indices
        return chain_copy

# A structure is a group of atoms organized in chains and residues
class Structure:
    def __init__ (self,
        atoms : List['Atom'] = [],
        residues : List['Residue'] = [],
        chains : List['Chain'] = [],
        ):
        self.atoms = []
        self.residues = []
        self.chains = []
        # Set references between instances
        for atom in atoms:
            self.set_new_atom(atom)
        for residue in residues:
            self.set_new_residue(residue)
        for chain in chains:
            self.set_new_chain(chain)
        # --- Set other internal variables ---
        # Set bonds between atoms
        self._bonds = None
        # --- Set other auxiliar variables ---
        # Trajectory atom sorter is a function used to sort coordinates in a trajectory file
        # This function is generated when sorting indices in the structure
        self.trajectory_atom_sorter = None

    def __repr__ (self):
        return '<Structure (' + str(len(self.atoms)) + ' atoms)>'

    # The bonds between atoms (read only)
    def get_bonds (self) -> List[ List[int] ]:
        # Return the stored value, if exists
        if self._bonds:
            return self._bonds
        # If not, we must calculate the bonds using vmd
        self._bonds = self.get_covalent_bonds()
        return self._bonds
    # Force specific bonds
    def set_bonds (self, bonds : List[ List[int] ]):
        self._bonds = bonds
    bonds = property(get_bonds, set_bonds, None, "The structure bonds (read only)")

    # Set a new atom in the structure
    def set_new_atom (self, atom : 'Atom'):
        atom._structure = self
        new_atom_index = len(self.atoms)
        self.atoms.append(atom)
        atom._index = new_atom_index

    # Set a new residue in the structure
    # WARNING: Atoms must be set already before setting residues
    def set_new_residue (self, residue : 'Residue'):
        residue._structure = self
        new_residue_index = len(self.residues)
        self.residues.append(residue)
        residue._index = new_residue_index
        # In case the residue has atom indices, set relational indices on each atom
        for atom_index in residue.atom_indices:
            atom = self.atoms[atom_index]
            atom._residue_index = new_residue_index

    # Purge residue from the structure and its chain
    # This can be done only when the residue has no atoms left in the structure
    # Renumerate all residue indices which have been offsetted as a result of the purge
    def purge_residue (self, residue : 'Residue'):
        # Check the residue can be purged
        if residue not in self.residues:
            raise ValueError(str(residue) + ' is not in the structure already')
        if len(residue.atom_indices) > 0:
            raise ValueError(str(residue) + ' is still having atoms and thus it cannot be purged')
        # Get the current index of the residue to be purged
        purged_index = residue.index
        # Residues and their atoms below this index are not to be modified
        # Residues and their atoms over this index must be renumerated
        for affected_residue in self.residues[purged_index+1:]:
            # Chaging the index automatically changes all residues atoms '_residue_index' values
            # Chaging the index automatically changes its corresponding index in residue chain '_residue_indices'
            affected_residue.index -= 1
        # Finally, remove the current residue from the list of residues in the structure
        del self.residues[purged_index]

    # Set a new chain in the structure
    # WARNING: Residues and atoms must be set already before setting chains
    def set_new_chain (self, chain : 'Chain'):
        chain._structure = self
        new_chain_index = len(self.chains)
        self.chains.append(chain)
        chain._index = new_chain_index
        # In case the chain has residue indices, set relational indices on each residue
        for residue_index in chain.residue_indices:
            residue = self.residues[residue_index]
            residue._chain_index = new_chain_index

    # Purge chain from the structure
    # This can be done only when the chain has no residues left in the structure
    # Renumerate all chain indices which have been offsetted as a result of the purge
    def purge_chain (self, chain : 'Chain'):
        # Check the chain can be purged
        if chain not in self.chains:
            raise ValueError('Chain ' + chain.name + ' is not in the structure already')
        if len(chain.residue_indices) > 0:
            raise ValueError('Chain ' + chain.name + ' is still having residues and thus it cannot be purged')
        # Get the current index of the chain to be purged
        purged_index = chain.index
        # Chains and their residues below this index are not to be modified
        # Chains and their residues over this index must be renumerated
        for affected_chain in self.chains[purged_index+1:]:
            # Chaging the index automatically changes all chain residues '_chain_index' values
            affected_chain.index -= 1
        # Finally, remove the current chain from the list of chains in the structure
        del self.chains[purged_index]

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
            parsed_atom = Atom(name=name, element=element, coords=coords)
            parsed_atoms.append(parsed_atom)
        # Parse residues
        for prody_residue in prody_residues:
            name = prody_residue.getResname()
            number = int(prody_residue.getResnum())
            icode = prody_residue.getIcode()
            parsed_residue = Residue(name=name, number=number, icode=icode)
            atom_indices = [ int(index) for index in prody_residue.getIndices() ]
            parsed_residue.atom_indices = atom_indices
            parsed_residues.append(parsed_residue)
        # Parse chains
        for prody_chain in prody_chains:
            name = prody_chain.getChid()
            parsed_chain = Chain(name=name)
            residue_indices = [ int(residue.getResindex()) for residue in prody_chain.iterResidues() ]
            parsed_chain.residue_indices = residue_indices
            parsed_chains.append(parsed_chain)
        return cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)

    # Set the structure from a pdb file
    @classmethod
    def from_pdb_file(cls, pdb_filename : str):
        if not os.path.exists(pdb_filename):
            raise SystemExit('File "' + pdb_filename + '" not found')
        if not get_format(pdb_filename) == 'pdb':
            raise SystemExit('"' + pdb_filename + '" is not a name for a pdb file')
        # Read the pdb file line by line and set the parsed atoms, residues and chains
        parsed_atoms = []
        parsed_residues = []
        parsed_chains = []
        atom_index = -1
        residue_index = -1
        with open(pdb_filename, 'r') as file:
            for line in file:
                # Parse atoms only
                start = line[0:6]
                is_atom = start == 'ATOM  ' or start == 'HETATM'
                if not is_atom:
                    continue
                # Mine all atom data
                atom_name = line[11:16].strip()
                residue_name = line[17:21].strip()
                chain = line[21:22]
                residue_number = int(line[22:26])
                icode = line[26:27]
                if icode == ' ':
                    icode = ''
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                element = line[76:78].strip()
                # Set the parsed atom, residue and chain
                parsed_atom = Atom(name=atom_name, element=element, coords=(x_coord, y_coord, z_coord))
                parsed_residue = Residue(name=residue_name, number=residue_number, icode=icode)
                parsed_chain = Chain(name=chain)
                # Add the parsed atom to the list and update the current atom index
                parsed_atoms.append(parsed_atom)
                atom_index += 1
                # Check if we are in the same chain/residue than before
                same_chain = parsed_chains and parsed_chains[-1] == parsed_chain
                same_residue = same_chain and parsed_residues and parsed_residue.same_inputs_as(parsed_residues[-1])
                # Update the residue atom indices
                # If the residue equals the last parsed residue then use the previous instead
                if same_residue:
                    parsed_residue = parsed_residues[-1]
                    parsed_residue.atom_indices.append(atom_index)
                    # If it is the same residue then it will be the same chain as well so we can proceed
                    continue
                # Otherwise, include the new residue in the list and update the current residue index
                parsed_residues.append(parsed_residue)
                residue_index += 1
                parsed_residue.atom_indices.append(atom_index)
                # If the chain equals the last parsed chain then use the previous instead
                if same_chain:
                    parsed_chain = parsed_chains[-1]
                    parsed_chain.residue_indices.append(residue_index)
                    continue
                # Otherwise, include the new chain in the list
                parsed_chains.append(parsed_chain)
                parsed_chain.residue_indices.append(residue_index)
        return cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)

    # Fix atom elements by gueesing them when missing
    # Set all elements with the first letter upper and the second (if any) lower
    # Also check if atom elements are coherent with atom names
    # If 'trust' is set as False then we impose elements according to what we can guess from the atom name
    # Return True if any element was modified or False if not
    def fix_atom_elements (self, trust : bool = True) -> bool:
        modified = False
        for atom in self.atoms:
            # Make sure elements have the first letter cap and the second letter not cap
            if atom.element:
                new_element = first_cap_only(atom.element)
                if atom.element != new_element:
                    atom.element = new_element
                    modified = True
                # Check the element to match what we would guess from the atom name
                # In case it does not just warn the user
                guess = guess_name_element(atom.name)
                if atom.element != guess:
                    print('WARNING: Suspicious element for atom ' + atom.name + ' -> ' + atom.element + " (shoudn't it be " + guess + "?)")
                    if not trust:
                        atom.element = guess
                        modified = True
            # If elements are missing then guess them from atom names
            else:
                atom.element = guess_name_element(atom.name)
                modified = True
        if modified:
            print('WARNING: Atom elements have been modified')
        return modified

    # Generate a pdb file with current structure
    def generate_pdb_file(self, pdb_filename : str):
        with open(pdb_filename, "w") as file:
            file.write('REMARK mdtoolbelt generated pdb file\n')
            for a, atom in enumerate(self.atoms):
                residue = atom.residue
                index = str((a+1) % 100000).rjust(5)
                name = ' ' + atom.name.ljust(3) if len(atom.name) < 4 else atom.name
                residue_name = residue.name.ljust(4) if residue else 'XXX'.ljust(4)
                chain = atom.chain
                chain_name = atom.chain.name.rjust(1) if chain else 'X'
                residue_number = str(residue.number).rjust(4) if residue else '0'.rjust(4)
                icode = residue.icode.rjust(1) if residue else ' '
                coords = atom.coords
                x_coord, y_coord, z_coord = [ "{:.3f}".format(coord).rjust(8) for coord in coords ]
                occupancy = '1.00' # Just a placeholder
                temp_factor = '0.00' # Just a placeholder
                element = atom.element.rjust(2)
                atom_line = ('ATOM  ' + index + ' ' + name + ' ' + residue_name
                    + chain_name + residue_number + icode + '   ' + x_coord + y_coord + z_coord
                    + '  ' + occupancy + '  ' + temp_factor + '          ' + element).ljust(80) + '\n'
                file.write(atom_line)

    # Get the structure equivalent prody topology
    def get_prody_topology (self):
        # Generate the prody topology
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        prody_topology = prody.parsePDB(pdb_filename)
        os.remove(pdb_filename)
        return prody_topology

     # Get the structure equivalent pytraj topology
    def get_pytraj_topology (self):
        # Generate a pdb file from the current structure to feed pytraj
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        pytraj_topology = pytraj.load_topology(filename = pdb_filename)
        os.remove(pdb_filename)
        return pytraj_topology

    # Select atoms from the structure thus generating an atom indices list
    # Different tools may be used to make the selection:
    # - vmd (default)
    # - prody
    # - pytraj
    def select (self, selection_string : str, syntax : str = 'vmd') -> Optional['Selection']:
        if syntax == 'vmd':
            # Generate a pdb for vmd to read it
            pdb_filename = '.structure.pdb'
            self.generate_pdb_file(pdb_filename)
            # Use vmd to find atom indices
            atom_indices = get_vmd_selection_atom_indices(pdb_filename, selection_string)
            os.remove(pdb_filename)
            if len(atom_indices) == 0:
                return None
            return Selection(atom_indices)
        if syntax == 'prody':
            prody_topology = self.get_prody_topology()
            prody_selection = prody_topology.select(selection_string)
            if not prody_selection:
                return None
            return Selection.from_prody(prody_selection)
        if syntax == 'pytraj':
            pytraj_topology = self.get_pytraj_topology()
            pytraj_selection = pytraj_topology[selection_string]
            atom_indices = [ atom.index for atom in prody_selection.atoms ]
            if len(atom_indices) == 0:
                return None
            return Selection(atom_indices)

        print('WARNING: Syntax ' + syntax + ' is not supported')
        return None

    # Set a function to make selections using atom indices
    def select_atom_indices (self, atom_indices : List[int]) -> 'Selection':
        # Check atom indices to be in the structure
        atom_count = len(self.atoms)
        for atom_index in atom_indices:
            if atom_index >= atom_count:
                raise SystemExit('Atom index ' + str(atom_index) + ' is out of range (' + str(atom_count) + ')')
        return Selection(atom_indices)

    # Set a function to make selections using residue indices
    def select_residue_indices (self, residue_indices : List[int]) -> 'Selection':
        atom_indices = sum([ self.residues[index].atom_indices for index in residue_indices ], [])
        return Selection(atom_indices)

    # Get a selection with all atoms
    def select_all (self) -> 'Selection':
        atom_count = len(self.atoms)
        return Selection(list(range(atom_count)))

    # Invert a selection
    def invert_selection (self, selection : 'Selection') -> 'Selection':
        return Selection([ atom_index for atom_index in range(len(self.atoms)) if atom_index not in selection.atom_indices ])
    
    # Given a selection, get a list of residue indices for residues implicated
    # Note that if a single atom from the residue is in the selection then the residue index is returned
    def get_selection_residue_indices (self, selection : 'Selection') -> List[int]:
        return list(set([ self.atoms[atom_index].residue_index for atom_index in selection.atom_indices ]))

    # Create a new structure from the current using a selection to filter atoms
    def filter (self, selection : 'Selection') -> 'Structure':
        if not selection:
            raise SystemExit('No selection was passed')
        new_atoms = []
        new_residues = []
        new_chains = []
        # Get the selected atoms
        # Generate dictionaries with new indexes as keys and previous indexes as values for atoms, residues and chains
        # This is done with this structure for the residues and chains further to find the new indices fast
        new_atom_indices = {}
        # Collect also original indices to related atom residues and chains
        original_atom_residue_indices = []
        for new_index, original_index in enumerate(selection.atom_indices):
            # Make a copy of the selected atom in order to not modify the original one
            original_atom = self.atoms[original_index]
            new_atom = Atom(
                name=original_atom.name,
                element=original_atom.element,
                coords=original_atom.coords,
            )
            new_atoms.append(new_atom)
            # Save old and new indices in order to reconfigure all indices later
            new_atom_indices[original_index] = new_index
            original_atom_residue_indices.append(original_atom.residue_index)
        # Find the selected residues
        selected_residue_indices = list(set(original_atom_residue_indices))
        # Repeat the original/new indices backup we did before
        new_residue_indices = {}
        original_residue_atom_indices = []
        original_residue_chain_indices = []
        for new_index, original_index in enumerate(selected_residue_indices):
            # Make a copy of the selected residue in order to not modify the original one
            original_residue = self.residues[original_index]
            new_residue = Residue(
                name=original_residue.name,
                number=original_residue.number,
                icode=original_residue.icode,
            )
            new_residues.append(new_residue)
            # Save old and new indices in order to reconfigure all indices later
            new_residue_indices[original_index] = new_index
            original_residue_atom_indices.append(original_residue.atom_indices)
            original_residue_chain_indices.append(original_residue.chain_index)
        # Find the selected chains
        selected_chain_indices = list(set(original_residue_chain_indices))
        # Repeat the original/new indices backup we did before
        new_chain_indices = {}
        original_chain_residue_indices = []
        for new_index, original_index in enumerate(selected_chain_indices):
            # Make a copy of the selected chain in order to not modify the original one
            original_chain = self.chains[original_index]
            new_chain = Chain(
                name=original_chain.name,
            )
            new_chains.append(new_chain)
            # Save old and new indices in order to reconfigure all indices later
            new_chain_indices[original_index] = new_index
            original_chain_residue_indices.append(original_chain.residue_indices)
        # Finally, reset indices in all instances
        for a, atom in enumerate(new_atoms):
            atom.residue_index = new_residue_indices[ original_atom_residue_indices[a] ]
        for r, residue in enumerate(new_residues):
            atom_indices = []
            for original_index in original_residue_atom_indices[r]:
                new_index = new_atom_indices.get(original_index, None)
                if new_index != None:
                    atom_indices.append(new_index)
            residue.atom_indices = atom_indices
            residue.chain_index = new_chain_indices[ original_residue_chain_indices[r] ]
        for c, chain in enumerate(new_chains):
            residue_indices = []
            for original_index in original_chain_residue_indices[c]:
                new_index = new_residue_indices.get(original_index, None)
                if new_index != None:
                    residue_indices.append(new_index)
            chain.residue_indices = residue_indices
        return Structure(atoms=new_atoms, residues=new_residues, chains=new_chains)

    # Set chains on demand
    # If no selection is passed then the whole structure will be affected
    # If no chain is passed then a "chain by fragment" logic will be applied
    def chainer (self, selection : Optional['Selection'] = None, letter : Optional[str] = None):
        # If there is no selection we consider all atoms
        if not selection:
            selection = self.select_all()
        # If a letter is specified then the logic is way simpler
        if letter:
            self.set_selection_chain_name(selection, letter)
            return
        # If a letter is not specified we run the "fragments" logic
        fragments = self.find_fragments(selection)
        for fragment in fragments:
            chain_name = self.get_next_available_chain_name()
            if not chain_name:
                raise SystemExit('ERROR: There are more chains than letters in the alphabet')
            fragment_selection = Selection(fragment)
            self.set_selection_chain_name(fragment_selection, chain_name)

    # This is an alternative system to find protein chains (anything else is chained as 'X')
    # This system does not depend on VMD
    # It totally overrides previous chains since it is expected to be used only when chains are missing
    def raw_protein_chainer (self):
        current_chain = self.get_next_available_chain_name()
        previous_alpha_carbon = None
        for residue in self.residues:
            alpha_carbon = next((atom for atom in residue.atoms if atom.name == 'CA'), None)
            if not alpha_carbon:
                residue.set_chain('X')
                continue
            # Connected aminoacids have their alpha carbons at a distance of around 3.8 Ångstroms
            residues_are_connected = previous_alpha_carbon and calculate_distance(previous_alpha_carbon, alpha_carbon) < 4
            if not residues_are_connected:
                current_chain = self.get_next_available_chain_name()
            residue.set_chain(current_chain)
            previous_alpha_carbon = alpha_carbon

    # Find fragments* in a selection of atoms
    # * A fragment is a group of colvalently bond atoms
    # A list of fragments is returned, where every fragment is a list of atom indices
    # All atoms are searched if no selection is provided
    def find_fragments (self, selection : Optional['Selection'] = None) -> List[ List[int] ]:
        # If there is no selection we consider all atoms
        if not selection:
            selection = self.select_all()
        # Find covalent bonds between atoms
        covalent_bonds = self.get_covalent_bonds(selection)
        atom_indexed_covalent_bonds = { atom_index: covalent_bonds[i] for i, atom_index in enumerate(selection.atom_indices) }
        # Group the connected atoms in "fragments"
        fragments = []
        while len(atom_indexed_covalent_bonds.keys()) > 0:
            start_atom_index, bonds = list(atom_indexed_covalent_bonds.items())[0]
            del atom_indexed_covalent_bonds[start_atom_index]
            fragment_atom_indices = [ start_atom_index ]
            while len(bonds) > 0:
                next_atom_index = bonds[0]
                next_bonds = atom_indexed_covalent_bonds.get(next_atom_index, None)
                # If this atom is out of the selection then skip it
                if not(next_bonds):
                    bonds.remove(next_atom_index)
                    continue
                next_new_bonds = [ bond for bond in next_bonds if bond not in fragment_atom_indices + bonds ]
                bonds += next_new_bonds
                fragment_atom_indices.append(next_atom_index)
                bonds.remove(next_atom_index)
                del atom_indexed_covalent_bonds[next_atom_index]
            fragments.append(fragment_atom_indices)
        return fragments

    # Given an atom selection, set the chain for all these atoms
    # Note that the chain is changed in every whole residue, no matter if only one atom was selected
    def set_selection_chain_name (self, selection : 'Selection', letter : str):
        # Find if the letter belongs to an already existing chain
        chain = self.get_chain_by_name(letter)
        # If the chain does not exist yet then create it
        if not chain:
            chain = Chain(name=letter)
            self.set_new_chain(chain)
        # Convert the selection to residues
        # WARNING: Remember to never use a set when handling residues or you may eliminate 'duplicated' residues
        residue_indices = self.get_selection_residue_indices(selection)
        # Set the chain index of every residue to the new chain
        for residue_index in residue_indices:
            # WARNING: Chain index has to be read every iteration since it may change. Do not save it
            residue = self.residues[residue_index]
            residue.set_chain_index(chain.index)

    # Get the next available chain name
    # Find alphabetically the first letter which is not yet used as a chain name
    # If all letters in the alphabet are used already then return None
    def get_next_available_chain_name (self) -> Optional[str]:
        current_chain_names = [ chain.name for chain in self.chains ]
        return next((name for name in available_caps if name not in current_chain_names), None)

    # Get a chain by its name
    def get_chain_by_name (self, name : str) -> 'Chain':
        return next((c for c in self.chains if c.name == name), None)

    # Get a summary of the structure
    def display_summary (self):
        print('Atoms: ' + str(len(self.atoms)))
        print('Residues: ' + str(len(self.residues)))
        print('Chains: ' + str(len(self.chains)))
        for chain in self.chains:
            print('Chain ' + chain.name + ' (' + str(len(chain.residue_indices)) + ' residues)')
            print(' -> ' + chain.get_sequence())

    # There may be chains which are equal in the structure (i.e. same chain name)
    # This means we have a duplicated/splitted chain
    # Repeated chains are usual and they are usually supported but with some problems
    # Also, repeated chains ususally come with repeated residues, which means more problems (see explanation below)
    
    # Check repeated chains
    # Rename repeated chains if the fix_chains argument is True
    # WARNING: The fix is possible only if there are less chains than the number of letters in the alphabet
    # Although there is no limitation in this code for chain names, setting long chain names is not compatible with pdb format
    # Return True if there were any repeats
    def check_repeated_chains (self, fix_chains : bool = False, display_summary : bool = False) -> bool:
        # Order chains according to their names
        # Save also those chains which have a previous duplicate
        name_chains = {}
        repeated_chains = []
        for chain in self.chains:
            chain_name = chain.name
            current_name_chains = name_chains.get(chain_name, None)
            if not current_name_chains:
                name_chains[chain_name] = [chain]
            else:
                name_chains[chain_name].append(chain)
                repeated_chains.append(chain)
        # Display the summary of repeated chains if requested
        if display_summary:
            if len(repeated_chains) > 0:
                print('WARNING: There are repeated chains:')
                for chain_name, chains in name_chains.items():
                    chains_count = len(chains)
                    if chains_count > 1:
                        print('- Chain ' + chain_name + ' has ' + str(chains_count) + ' repeats' )
        # Rename repeated chains if requested
        if len(repeated_chains) > 0 and fix_chains:
            if len(self.chains) > 26:
                # for chain in self.chains:
                #     print(str(chain) + ': ' + str(chain.atom_indices[0]) + ' to ' + str(chain.atom_indices[-1]))
                raise ValueError('There are more chains than letters in the alphabet')
            current_letters = list(name_chains.keys())
            for repeated_chain in repeated_chains:
                last_chain_letter = repeated_chain.name
                while last_chain_letter in current_letters:
                    last_chain_letter = get_next_letter(last_chain_letter)
                repeated_chain.name = last_chain_letter
                current_letters.append(last_chain_letter)
        # Fix repeated chains if requested
        return len(repeated_chains) > 0

    # There may be residues which are equal in the structure (i.e. same chain, name, number and icode)
    # In case 2 residues in the structure are equal we must check distance between their atoms
    # If atoms are far it means they are different residues with the same notation (duplicated residues)
    # If atoms are close it means they are indeed the same residue (splitted residue)

    # Splitted residues are found in some pdbs and they are supported by some tools
    # These tools consider all atoms with the same 'record' as the same residue
    # However, there are other tools which would consider the splitted residue as two different resdiues
    # This causes inconsistency along different tools besides a long list of problems
    # This situation has no easy fix since we can not change the order in atoms because trajectory data depends on it
    # DANI: Habrá que hacer una función para reordenar átomos en la estructura a la vez que en una trayectoria

    # Duplicated residues are usual and they are usually supported but with some problems
    # For example, pytraj analysis outputs use to sort results by residues and each residue is tagged
    # If there are duplicated residues with the same tag it may be not possible to know which result belongs to each residue
    # Another example are NGL selections once in the web client
    # If you select residue ':A and 1' and there are multiple residues 1 in chain A all of them will be displayed

    # Check residues to search for duplicated and splitted residues
    # Renumerate repeated residues if the fix_residues argument is True
    # Return True if there were any repeats
    def check_repeated_residues (self, fix_residues : bool = False, display_summary : bool = False) -> bool:
        # Track if residues have to be changed or not
        modified = False
        # Group all residues in the structure according to their chain, number and icode
        grouped_residues = {}
        # Check repeated residues which are one after the other
        # Note that these residues MUST have a different name
        # Otherwise they would have not been considered different residues
        # For these rare cases we use icodes to solve the problem
        non_icoded_residues = []
        last_residue = None
        for residue in self.residues:
            # Check residue to be equal than the previous residue
            if residue == last_residue:
                non_icoded_residues.append(residue)
                last_residue = residue
                continue
            last_residue = residue
            # Add residue to the groups of repeated residues
            current_residue_repeats = grouped_residues.get(residue, None)
            if not current_residue_repeats:
                grouped_residues[residue] = [ residue ]
            else:
                current_residue_repeats.append(residue)
        # In case we have non icoded residues
        if len(non_icoded_residues) > 0:
            if display_summary:
                print('There are non-icoded residues (' + str(len(non_icoded_residues)) + ')')
            # Set new icodes for non icoded residues
            if fix_residues:
                print('    Non icoded residues will recieve an icode')
                for residue in non_icoded_residues:
                    repeated_residues_group = grouped_residues[residue]
                    current_icodes = [ residue.icode for residue in repeated_residues_group if residue.icode ]
                    next_icode = next((cap for cap in available_caps if cap not in current_icodes), None)
                    if not next_icode:
                        raise ValueError('There are no more icodes available')
                    # print('Setting icode ' + next_icode + ' to residue ' + str(residue))
                    residue.icode = next_icode
                modified = True
        # Grouped resdiues with more than 1 result are considered as repeated
        repeated_residues = [ residues for residues in grouped_residues.values() if len(residues) > 1 ]
        if len(repeated_residues) == 0:
            return modified
        # In case we have repeated residues...
        if display_summary:
            print('WARNING: There are repeated residues (' + str(len(repeated_residues)) + ')')
            print('    e.g. ' + str(repeated_residues[0][0]))
        # Now for each repeated residue, find out which are splitted and which are duplicated
        covalent_bonds = self.get_covalent_bonds()
        overall_splitted_residues = []
        overall_duplicated_residues = []
        for residues in repeated_residues:
            # Iterate over repeated residues and check if residues are covalently bonded
            # If any pair of residues are bonded add them both to the splitted residues list
            # At the end, all non-splitted residues will be considered duplicated residues
            # WARNING: Using a set here is not possible since repeated residues have the same hash
            # WARNING: Also comparing residues themselves is not advisable, so we use indices at this point
            splitted_residue_indices = set()
            for residue, other_residues in otherwise(residues):
                if residue.index in splitted_residue_indices:
                    continue
                # Get atom indices for all atoms connected to the current residue
                residue_bonds = sum([ covalent_bonds[index] for index in residue.atom_indices ], [])
                for other_residue in other_residues:
                    # Get all atom indices for each other residue and collate with the current residue bonds
                    if any( index in residue_bonds for index in other_residue.atom_indices ):
                        splitted_residue_indices.add(residue.index)
                        splitted_residue_indices.add(other_residue.index)
            # Finally obtain the splitted residues from its indices
            splitted_residues = [ self.residues[index] for index in splitted_residue_indices ]
            # Repeated residues which are not splitted are thus duplicated
            duplicated_residues = [ residue for residue in residues if residue.index not in splitted_residue_indices ]
            # Update the overall lists
            if len(splitted_residues) > 0:
                overall_splitted_residues.append(splitted_residues)
            if len(duplicated_residues) > 0:
                overall_duplicated_residues += duplicated_residues
        # In case we have splitted residues
        if len(overall_splitted_residues) > 0:
            if display_summary:
                print('    There are splitted residues (' + str(len(overall_splitted_residues)) + ')')
            # Fix splitted residues
            if fix_residues:
                print('        Splitted residues will be merged')
                # Set a function to sort atoms and residues by index
                def by_index (v):
                    return v._index
                # Merge splitted residues
                # WARNING: Merging residues without sorting atoms is possible, but this would be lost after exporting to pdb
                for splitted_residues in overall_splitted_residues:
                    # The first residue (i.e. the residue with the lower index) will be the one which remains
                    # It will absorb other residue atom indices
                    splitted_residues.sort(key=by_index) # Residues are not sorted by default, this is mandatory
                    first_residue = splitted_residues[0]
                    other_residues = splitted_residues[1:]
                    for residue in other_residues:
                        for atom in residue.atoms:
                            atom.residue = first_residue
                        # Then these other residues are eliminated
                        self.purge_residue(residue)
                print('        Atoms will be sorted to be together by residues')
                print('        NEVER FORGET: This will break any associated trajectory if coordinates are not sorted as well')
                # Sort atoms to group residue atoms together
                # Note that each atom index must be updated
                new_atom_indices = sum([ residue.atom_indices for residue in self.residues ], [])
                for i, new_atom_index in enumerate(new_atom_indices):
                    atom = self.atoms[new_atom_index]
                    atom._index = i
                self.atoms.sort(key=by_index)
                # Also residue 'atom_indices' must be updated
                for residue in self.residues:
                    residue._atom_indices = [ new_atom_indices.index(atom_index) for atom_index in residue._atom_indices ]
                # Prepare the trajectory atom sorter which must be returned
                # Include atom indices already so the user has to provide only the structure and trajectory filenames
                def trajectory_atom_sorter (
                    input_structure_filename : str,
                    input_trajectory_filename : str,
                    output_trajectory_filename : str
                ):
                    sort_trajectory_atoms(
                        input_structure_filename,
                        input_trajectory_filename,
                        output_trajectory_filename,
                        new_atom_indices
                    )
                self.trajectory_atom_sorter = trajectory_atom_sorter
                modified = True
         # In case we have duplicated residues
        if len(overall_duplicated_residues) > 0:
            if display_summary:
                print('    There are duplicated residues (' + str(len(overall_duplicated_residues)) + ')')
            # Renumerate duplicated residues if requested
            if fix_residues:
                print('        Duplicated residues will be renumerated')
                # Get the next available number in the residue chain
                for duplicated_residue in overall_duplicated_residues:
                    maximum_chain_number = max([ residue.number for residue in duplicated_residue.chain.residues ])
                    duplicated_residue.number = maximum_chain_number + 1
                modified = True
        return modified

    # Atoms with identical chain, residue and name are considered repeated atoms
    # DANI: No recuerdo los problemas que daba tener átomos repetidos

    # Check atoms to search for repeated atoms
    # Rename repeated atoms if the fix_atoms argument is True
    # Return True if there were any repeats
    def check_repeated_atoms (self, fix_atoms : bool = False, display_summary : bool = False) -> bool:
        # Set two trackers for display
        repeated_atoms_count = 0
        example = None
        for residue in self.residues:
            # Iterate over the residue atoms counting how many repeated names we have
            current_names = []
            for atom in residue.atoms:
                # We check if the atom name already exists. If not, go to the next atom
                name = atom.name
                if name not in current_names:
                    current_names.append(name)
                    continue
                # When atom is repeated
                repeated_atoms_count += 1
                # If the fix was not requested we stop here
                if not fix_atoms:
                    continue
                # We set the name of the atom as the element letter + the count of this element
                # If element is missing take the first character of the atom name
                initial = atom.element
                if not initial or initial == ' ':
                    initial = name[0]
                number = 1
                new_name = initial + str(number)
                while new_name in current_names:
                    number += 1
                    new_name = initial + str(number)
                # Save an example for the logs if there is none yet
                if not example:
                    example = atom.name + ' renamed as ' + new_name + ' in residue ' + str(residue)
                atom.name = new_name
                current_names.append(new_name)
        # Display the summary of repeated atoms if requested
        if display_summary:
            if repeated_atoms_count > 0:
                print('WARNING: There are repeated atoms (' + str(repeated_atoms_count) + ')')
                print('    e.g. ' + example)
        return repeated_atoms_count > 0

    # Get all atomic covalent (strong) bonds
    # Bonds are defined as a list of atom indices for each atom in the structure
    # Rely on VMD logic to do so
    def get_covalent_bonds (self, selection : Optional['Selection'] = None) -> List[ List[int] ]:
        # Generate a pdb strucutre to feed vmd
        auxiliar_pdb_filename = '.structure.pdb'
        self.generate_pdb_file(auxiliar_pdb_filename)
        # Get covalent bonds between both residue atoms
        covalent_bonds = get_covalent_bonds(auxiliar_pdb_filename, selection)
        # Remove the auxiliar pdb file
        os.remove(auxiliar_pdb_filename)
        return covalent_bonds

    # Make a copy of the current structure
    def copy (self) -> 'Structure':
        atom_copies = [ atom.copy() for atom in self.atoms ]
        residue_copies = [ residue.copy() for residue in self.residues ]
        chain_copies = [ chain.copy() for chain in self.chains ]
        return Structure(atom_copies, residue_copies, chain_copies)

    # Merge currnet structure with another structure
    # DANI: No lo he testeado en profundidad
    def merge (self, other : 'Structure') -> 'Structure':
        # Copy self atoms, residues and chains
        self_atom_copies = [ atom.copy() for atom in self.atoms ]
        self_residue_copies = [ residue.copy() for residue in self.residues ]
        self_chain_copies = [ chain.copy() for chain in self.chains ]
        # Copy other atoms, residues and chains
        other_atom_copies = [ atom.copy() for atom in other.atoms ]
        other_residue_copies = [ residue.copy() for residue in other.residues ]
        other_chain_copies = [ chain.copy() for chain in other.chains ]
        # Adapt indices in other atoms, residues and chains
        atom_index_offset = len(self_atom_copies)
        residue_index_offset = len(self_residue_copies)
        chain_index_offset = len(self_chain_copies)
        for atom in other_atom_copies:
            atom._index += atom_index_offset
            atom._residue_index += residue_index_offset
        for residue in other_residue_copies:
            residue._index += residue_index_offset
            residue._atom_indices = [ i + atom_index_offset for i in residue._atom_indices ]
            residue._chain_index += chain_index_offset
        for chain in other_chain_copies:
            chain._index += chain_index_offset
            chain._residue_indices = [ i + residue_index_offset for i in chain._residue_indices ]
        # Merge self with other atoms, residues and chains and build the new structure
        merged_atoms = self_atom_copies + other_atom_copies
        merged_residues = self_residue_copies + other_residue_copies
        merged_chains = self_chain_copies + other_chain_copies
        return Structure(merged_atoms, merged_residues, merged_chains)
            

### Related functions ###

# Calculate the distance between two atoms
def calculate_distance (atom_1 : Atom, atom_2 : Atom) -> float:
    squared_distances_sum = 0
    for i in { 'x': 0, 'y': 1, 'z': 2 }.values():
        squared_distances_sum += (atom_1.coords[i] - atom_2.coords[i])**2
    return math.sqrt(squared_distances_sum)

### Auxiliar functions ###

# Set a function to get the next letter from an input letter in alphabetic order
def get_next_letter (letter : str) -> str:
    if letter == 'z' or letter == 'Z':
        raise ValueError("Limit of chain letters has been reached")
    next_letter = chr(ord(letter) + 1)
    return next_letter

# Given a string with 1 or 2 characters, return a new string with the first letter cap and the second letter not cap (if any)
def first_cap_only (name : str) -> str:
    if len(name) == 1:
        return name.upper()
    first_character = name[0].upper()
    second_character = name[1].lower()
    return first_character + second_character

# Guess an atom element from its name
# DANI: Cuidado con añadir el calcio (Ca) a la lista, porque te cambia todo los carbonos alpha
supported_elements = [ 'C', 'N', 'O', 'H', 'P', 'S', 'K', 'F', 'Cl', 'Na', 'Zn', 'Mg', 'Fe', 'Br', 'Mn' ]
def guess_name_element (name: str) -> str:
    length = len(name)
    next_character = None
    for i, character in enumerate(name):
        # Get the next character, since element may be formed by 2 letters
        if i < length - 1:
            next_character = name[i+1]
            # If next character is not a string then ignore it
            if not next_character.isalpha():
                next_character = None
        # Try to get all possible matches between the characters and the supported atoms
        # First letter is always caps
        character = character.upper()
        # First try to match both letters together
        if next_character:
            # Start with the second letter in caps
            next_character = next_character.upper()
            both = character + next_character
            if both in supported_elements:
                return both
            # Continue with the second letter in lowers
            next_character = next_character.lower()
            both = character + next_character
            if both in supported_elements:
                return both
        # Finally, try with the first character alone
        if character in supported_elements:
            return character
    raise SystemExit(
        "ERROR: Not recognized element in '" + name + "'")

# Set a special iteration system
# Return one value of the array and a new array with all other values for each value
def otherwise (values : list) -> Generator[tuple, None, None]:
    for v, value in enumerate(values):
        others = values[0:v] + values[v+1:]
        yield value, others