"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a single neutral atom molecule is a metal atom based on its SMILES string.
    Metal atoms are identified from the periodic table, and they have specific atomic numbers.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a neutral metal atom, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for single atom (possibly isotopic)
    if mol.GetNumAtoms() != 1:
        return False, "Molecule does not consist of a single atom"

    # Get the single atom
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    isotope = atom.GetIsotope()  # Retrieve isotope information
    formal_charge = atom.GetFormalCharge()
    
    # Ensure the atom is neutral
    if formal_charge != 0:
        return False, f"Atom has a non-zero formal charge: {formal_charge}"

    # List of metal atomic numbers based on a comprehensive periodic table inclusion
    metal_atomic_numbers = {
        3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58,
        59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
        77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
        97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
        112, 113, 114, 115, 116, 117
    }
    
    # Check if atomic number is classified as a metal and the atom is neutral
    if atomic_num in metal_atomic_numbers:
        return True, f"Neutral atom with atomic number {atomic_num}{(' with isotope ' + str(isotope)) if isotope else ''} is a metal atom"
    else:
        return False, f"Atom with atomic number {atomic_num}{(' with isotope ' + str(isotope)) if isotope else ''} is not a metal atom"