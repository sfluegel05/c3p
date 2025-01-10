"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a molecule represented by a SMILES string is a metal atom.
    
    A metal atom is defined as an atom of an element that exhibits typical 
    metallic properties, being typically shiny, with high electrical and thermal 
    conductivity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the SMILES represents a metal atom, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not represent a single metal atom"

    # Get the atomic information
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    symbol = atom.GetSymbol()

    # Strip any isotope information by getting the symbol name without numeric components
    pure_symbol = ''.join([i for i in symbol if not i.isdigit()])

    # List of atomic numbers for metal elements based on metallic properties
    metal_atomic_nums = [
        3,  4,  11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55,
        56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
        73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90,
        91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
        105, 106, 107, 108, 110, 111, 112
    ]

    # Check if the atomic number corresponds to a metal
    if atomic_num in metal_atomic_nums:
        return True, f"The atom is a metal with symbol '{pure_symbol}' and atomic number {atomic_num}"

    return False, f"The atom with symbol '{pure_symbol}' and atomic number {atomic_num} is not a metal"

# Example usage:
# print(is_metal_atom("[23Na]"))
# print(is_metal_atom("[6C]"))  # Non-metal example