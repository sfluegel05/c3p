"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is defined as an atom of an element that exhibits typical metallic properties.

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

    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Molecule must consist of a single atom"

    # Get the atomic number of the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_number = atom.GetAtomicNum()

    # Check if the atom is neutral (charge 0)
    if atom.GetFormalCharge() != 0:
        return False, "Atom must be neutral (charge 0)"

    # List of atomic numbers corresponding to metals (including heavier elements)
    metal_atomic_numbers = {
        3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
        55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
        84,  # Polonium
        87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109  # Meitnerium
    }

    # Check if the atomic number is in the list of metals
    if atomic_number in metal_atomic_numbers:
        return True, f"Atom with atomic number {atomic_number} is a metal"
    else:
        return False, f"Atom with atomic number {atomic_number} is not a metal"