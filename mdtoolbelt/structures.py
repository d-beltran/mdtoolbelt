# Main handler of the toolbelt
import os
from bisect import bisect
from typing import Optional, Tuple, List, Union
Coords = Tuple[float, float, float]

import prody

from .selections import Selection
from .vmd_spells import get_vmd_selection_atom_indices

# ------------------------------------------------------------------------------------

# Functions for type handling

# Get the residue and the residue index out of both the residue or the residue index
def parse_residue (input_residue : Union['Residue', int]) -> Tuple['Residue', int]:
    if type(input_residue) == int:
        residue_index = input_residue
        residue = residue.structure.residues[residue_index]
    elif type(input_residue) == 'Residue':
        residue = input_residue
        residue_index = residue.index
    else:
        raise ValueError('Unknow type when expecting Residue or int')
    return residue, residue_index

# Get the chain and the chain index out of both the chain or the chain index
def parse_chain (input_chain : Union['Chain', int]) -> Tuple['Chain', int]:
    if type(input_chain) == int:
        chain_index = input_chain
        chain = chain.structure.chains[chain_index]
    elif type(input_chain) == 'Chain':
        chain = input_chain
        chain_index = chain.index
    else:
        raise ValueError('Unknow type when expecting Chain or int')
    return chain, chain_index

# ------------------------------------------------------------------------------------

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
        self.structure = None
        self.index = None
        self.residue = None
        self.chain = None
        self.residue_index = None
        self.chain_index = None

    def __repr__ (self):
        return '<Atom ' + self.name + '>'


# A residue
class Residue:
    def __init__ (self,
        name : Optional[str] = None,
        number : Optional[int] = None,
        icode : Optional[str] = None,
        atom_indices : List[int] = [],
        ):
        self.name = name
        self.number = number
        self.icode = icode
        self.atom_indices = atom_indices
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self.structure = None
        self.index = None
        self.atoms = []
        self.chain = None
        self.chain_index = None

    def __repr__ (self):
        return '<Residue ' + self.name + str(self.number) + (self.icode if self.icode else '') + '>'

    # Change the chain of the residue
    def change_chain (self, new_chain : Union['Chain', int]):
        # Get the residue current chain and remove the current residue from it
        current_chain = self.chain
        current_chain.remove_residue(self)
        # Get the new chain to be set
        new_chain, new_chain_index = parse_chain(new_chain)
        new_chain.add_residue(self)

# A chain
class Chain:
    def __init__ (self,
        name : Optional[str] = None,
        residue_indices : List[int] = [],
        ):
        self.name = name
        self.residue_indices = residue_indices
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self.structure = None
        self.index = None
        self.residues = []
        self.atoms = []
        self.atom_indices = []

    def __repr__ (self):
        return '<Chain ' + self.name + '>'

    # Add a residue to the chain
    def add_residue (self, input_residue : Union['Residue', int]):
        residue, residue_index = parse_residue(input_residue)
        sorted_residue_index = bisect(self.residue_indices, residue_index)
        self.residue_indices.insert(sorted_residue_index, residue_index)
        self.residues.insert(sorted_residue_index, residue)
        for atom_index in residue.atom_indices:
            sorted_atom_index = bisect(self.atom_indices, atom_index)
            self.atom_indices.insert(sorted_atom_index, atom_index)
            atom = self.structure.atoms[atom_index]
            self.atoms.insert(sorted_atom_index, atom)
            # Update atoms themselves
            atom.chain = self
            atom.chain_index = self.index
        # Update the residue itself
        residue.chain = self
        residue.chain_index = self.index

    # Remove a residue from the chain
    def remove_residue (self, input_residue : Union['Residue', int]):
        residue, residue_index = parse_residue(input_residue)
        self.residue_indices.remove(residue_index) # This index MUST be in the list
        self.residues.pop(residue_index)
        for atom_index in residue.atom_indices:
            self.atom_indices.remove(atom_index) # This index MUST be in the list
            self.atoms.pop(atom_index)
            # Update atoms themselves
            atom = self.structure.atoms[atom_index]
            atom.chain = None
            atom.chain_index = None
        # Update the residue itself
        residue.chain = None
        residue.chain_index = None

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
        atom.structure = self
        new_atom_index = len(self.atoms)
        self.atoms.append(atom)
        atom.index = new_atom_index

    # Set a new residue in the structure
    # WARNING: Atoms must be set already before setting residues
    def set_new_residue (self, residue : 'Residue'):
        residue.structure = self
        new_residue_index = len(self.residues)
        self.residues.append(residue)
        residue.index = new_residue_index
        for atom_index in residue.atom_indices:
            atom = self.atoms[atom_index]
            residue.atoms.append(atom)
            atom.residue_index = new_residue_index
            atom.residue = residue

    # Set a new chain in the structure
    # WARNING: Residues and atoms must be set already before setting chains
    def set_new_chain (self, chain : 'Chain'):
        chain.structure = self
        new_chain_index = len(self.chains)
        self.chains.append(chain)
        chain.index = new_chain_index
        # Set chain residues and update those residues about the chain they belong to
        for residue_index in chain.residue_indices:
            residue = self.residues[residue_index]
            chain.residues.append(residue)
            residue.chain_index = new_chain_index
            residue.chain = chain
            # Set chain atoms and update those atoms about the chain they belong to
            for atom_index in residue.atom_indices:
                atom = self.atoms[atom_index]
                chain.atoms.append(atom)
                atom.chain_index = new_chain_index
                atom.chain = chain

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
            atom_indices = [ int(index) for index in prody_residue.getIndices() ]
            parsed_residue = Residue(name=name, number=number, icode=icode, atom_indices=atom_indices)
            parsed_residues.append(parsed_residue)
        # Parse chains
        for prody_chain in prody_chains:
            name = prody_chain.getChid()
            residue_indices = [ int(residue.getResindex()) for residue in prody_chain.iterResidues() ]
            parsed_chain = Chain(name=name, residue_indices=residue_indices)
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
        original_atom_chain_indices = []
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
            original_atom_chain_indices.append(original_atom.chain_index)
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
        original_chain_atom_indices = []
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
            original_chain_atom_indices.append(original_chain.atom_indices)
            original_chain_residue_indices.append(original_chain.residue_indices)
        # Finally, reset indices in all instances
        for a, atom in enumerate(new_atoms):
            atom.residue_index = new_residue_indices[ original_atom_residue_indices[a] ]
            atom.chain_index = new_chain_indices[ original_atom_chain_indices[a] ]
        for r, residue in enumerate(new_residues):
            residue.atom_indices = []
            for original_index in original_residue_atom_indices[r]:
                new_index = new_atom_indices.get(original_index, None)
                if new_index != None:
                    residue.atom_indices.append(new_index)
            residue.chain_index = new_chain_indices[ original_residue_chain_indices[r] ]
        for c, chain in enumerate(new_chains):
            chain.atom_indices = []
            for original_index in original_chain_atom_indices[c]:
                new_index = new_atom_indices.get(original_index, None)
                if new_index != None:
                    chain.atom_indices.append(new_index)
            chain.residue_indices = []
            for original_index in original_chain_residue_indices[c]:
                new_index = new_residue_indices.get(original_index, None)
                if new_index != None:
                    chain.residue_indices.append(new_index)
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
        #     residue.change_chain(new_chain)

        pass

    # Get a chain by its name
    def get_chain_by_name (self, name : str) -> 'Chain':
        return next((c for c in self.chains if c.name == name), None)


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