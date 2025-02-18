"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

# List of atomic numbers for metallic elements (including alkali, alkaline earth, transition, 
# post-transition metals, lanthanides, actinides, and metalloids like polonium)
METAL_ATOMIC_NUMBERS = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55,
    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106,
    107, 108, 109, 110, 111, 112, 113, 114, 115, 116
}

def is_metal_atom(smiles: str):
    """
    Determines if a SMILES string represents a single metal atom in elemental state.
    Checks for atomic structure, metallic element, neutral charge, and no radicals.

    Args:
        smiles (str): Input SMILES string

    Returns:
        bool: True if metal atom, False otherwise
        str: Reason for decision
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # Check SMILES validity
    if not mol:
        return False, "Invalid SMILES"
    
    # Must be exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    atom = mol.GetAtomWithIdx(0)
    
    # Check if element is in metal list
    if atom.GetAtomicNum() not in METAL_ATOMIC_NUMBERS:
        return False, "Non-metallic element"
    
    # Validate neutral charge (elemental form)
    if atom.GetFormalCharge() != 0:
        return False, "Charged species"
    
    # Check for radicals (unlikely but possible)
    if atom.GetNumRadicalElectrons() != 0:
        return False, "Radical electrons present"
    
    return True, "Single metallic atom in elemental state"