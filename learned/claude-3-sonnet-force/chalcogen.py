"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHEBI:33373 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    A chalcogen is any p-block element belonging to the group 16 family of the periodic table.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the atomic numbers of all atoms in the molecule
    atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    
    # Check if the molecule contains any chalcogen elements
    chalcogens = [8, 16, 34, 52, 84, 118]  # Atomic numbers of chalcogen elements
    if any(atomic_number in chalcogens for atomic_number in atomic_numbers):
        return True, "The molecule contains one or more chalcogen elements"
    else:
        return False, "The molecule does not contain any chalcogen elements"