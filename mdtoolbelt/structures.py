# Main handler of the toolbelt
import os
import math
from bisect import bisect
from typing import Optional, Tuple, List
Coords = Tuple[float, float, float]

import prody

from .selections import Selection
from .vmd_spells import get_vmd_selection_atom_indices

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

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self):
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self):
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
    def get_residue (self) -> 'Residue':
        # If there is not strucutre yet it means the atom is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        if not self.structure:
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
    def get_chain_index (self) -> int:
        return self.residue.chain_index
    chain_index = property(get_chain_index, None, None, "The atom chain index according to parent structure chains (read only)")

    # The atom chain (read only)
    # In order to change the chain it must be changed in the atom residue
    def get_chain (self) -> 'Chain':
        # Get the chain in the structure according to the chain index
        return self.structure.chains[self.chain_index]
    chain = property(get_chain, None, None, "The atom chain (read only)")



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
        return self.name == other.name and self.number == other.number and self.icode == other.icode

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self):
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self):
        return self._index
    index = property(get_index, None, None, "The residue index according to parent structure residues (read only)")

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
        atom._residue_index = self.index

    # Remove an atom from the residue
    def remove_atom (self, current_atom : 'Atom'):
        # Remove the current atom index from the atom indices list
        self.atom_indices.remove(current_atom.index) # This index MUST be in the list
        # Update the atom internal index
        atom._residue_index = None

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
        current_chain = self.chain
        current_chain.remove_residue(self)
        new_chain = self.structure.chains[new_chain_index]
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
    def set_chain (self, new_chain : 'Chain'):
        # Find the new chain index and set it as the residue chain index
        # Note that the chain must be set in the structure already
        new_chain_index = new_chain.index
        if new_chain_index == None:
            raise ValueError('Chain ' + str(new_chain) + ' is not set in the structure')
        self.set_chain_index(new_chain_index)
    chain = property(get_chain, set_chain, None, "The residue chain")
    

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
        self.residue_indices = []

    def __repr__ (self):
        return '<Chain ' + self.name + '>'

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self):
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self):
        return self._index
    index = property(get_index, None, None, "The residue index according to parent structure residues (read only)")

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
        residue.chain_index = self.index

    # Remove a residue from the chain
    def remove_residue (self, residue : 'Residue'):
        self.residue_indices.remove(residue.index) # This index MUST be in the list
        # Update the residue internal chain index
        residue.chain_index = None

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

    def __repr__ (self):
        return '<Structure (' + str(len(self.atoms)) + ' atoms)>'

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
    # Use ProDy to do so
    @classmethod
    def from_pdb_file(cls, pdb_filename : str):
        if not os.path.exists(pdb_filename):
            raise SystemExit('File "' + pdb_filename + '" not found')
        prody_topology = prody.parsePDB(pdb_filename)
        return cls.from_prody(prody_topology)

    # Fix atom elements by gueesing them when missing
    # Set all elements with the first letter upper and the second (if any) lower
    def fix_atom_elements (self):
        for atom in self.atoms:
            # Make sure elements have the first letter cap and the second letter not cap
            if atom.element:
                atom.element = first_cap_only(atom.element)
            # If elements are missing then guess them from atom names
            else:
                atom.element = guess_name_element(atom.name)

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
                    + '  ' + occupancy + '  ' + temp_factor + '           ' + element).ljust(80) + '\n'
                file.write(atom_line)

    # Get the structure equivalent prody topology
    def get_prody_topology (self):
        # Generate the prody topology
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        prody_topology = prody.parsePDB(pdb_filename)
        os.remove(pdb_filename)
        return prody_topology

    # Select atoms from the structure thus generating an atom indices list
    # Different tools may be used to make the selection:
    # - vmd (default)
    # - prody
    def select (self, selection_string : str, syntax : str = 'vmd') -> Optional['Selection']:
        if syntax == 'vmd':
            # Generate a pdb for vmd to read it
            pdb_filename = '.structure.pdb'
            self.generate_pdb_file(pdb_filename)
            # Use vmd to find atom indices
            atom_indices = get_vmd_selection_atom_indices(pdb_filename, selection_string)
            os.remove(pdb_filename)
            if len(atom_indices) == 0:
                print('WARNING: Empty selection')
                return None
            return Selection(atom_indices)
        if syntax == 'prody':
            prody_topology = self.get_prody_topology()
            prody_selection = prody_topology.select(selection_string)
            if not prody_selection:
                print('WARNING: Empty selection')
                return None
            return Selection.from_prody(prody_selection)
        print('WARNING: Syntax ' + syntax + ' is not supported')
        return None
    
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

        # DANI: Aquí falta encontrar los strong bonds, para lo cual hay que extraer la lógica del workflow

        # DANI: Si la cadena no existe habrá que crearla
        #     new_chain = Chain(name=new_chain.name, residue_indices=[self.index])
        #     self.structure.set_new_chain(new_chain) # This function updates atoms and this residue already

        # DANI: Cuando el cambio esté claro para cada residuo:
        #     residue.chain = new_chain

        pass

    # Get a chain by its name
    def get_chain_by_name (self, name : str) -> 'Chain':
        return next((c for c in self.chains if c.name == name), None)

### Related functions ###

# Calculate the distance between two atoms
coordinate_indices = { 'x': 0, 'y': 1, 'z': 2 }
def caluclate_distance (atom_1 : Atom, atom_2 : Atom) -> float:
    squared_distances_sum = 0
    for i in coordinate_indices.values():
        squared_distances_sum += (atom_1.coords[i] - atom_2.coords[i])**2
    return math.sqrt(squared_distances_sum)


### Auxiliar functions ###

# Given a string with 1 or 2 characters, return a new string with the first letter cap and the second letter not cap (if any)
def first_cap_only (name : str) -> str:
    if len(name) == 1:
        return name.upper()
    first_character = name[0].upper()
    second_character = name[1].lower()
    return first_character + second_character

# Guess an atom element from its name
supported_elements = [ 'C', 'N', 'O', 'H', 'P', 'S', 'K', 'F', 'Cl', 'Na', 'Zn', 'Mg', 'Fe' ]
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