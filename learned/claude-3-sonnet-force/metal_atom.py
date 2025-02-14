"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: CHEBI:33441 metal atom
"""
from rdkit import Chem

METAL_ATOMIC_NUMBERS = set(range(21, 104))  # Atomic numbers of metals
METAL_ATOMIC_NUMBERS.update([112, 114])  # Add heavier metals Cd and In

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is an atom of an element that exhibits typical metallic properties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has only one atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"

    # Get atomic number of the atom
    atom = mol.GetAtoms()[0]
    atomic_num = atom.GetAtomicNum()

    # Check if atomic number corresponds to a metal
    if atomic_num in METAL_ATOMIC_NUMBERS:
        return True, f"Atomic number {atomic_num} corresponds to a metal atom"
    else:
        return False, f"Atomic number {atomic_num} does not correspond to a metal atom"


# Examples
print(is_metal_atom("[Y]"))  # True, 'Atomic number 39 corresponds to a metal atom'
print(is_metal_atom("[Na]"))  # True, 'Atomic number 11 corresponds to a metal atom'
print(is_metal_atom("[He]"))  # False, 'Atomic number 2 does not correspond to a metal atom'
print(is_metal_atom("[H2O]"))  # False, 'Not a single atom'
print(is_metal_atom("invalid_smiles"))  # False, 'Invalid SMILES string'