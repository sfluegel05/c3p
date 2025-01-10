"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a SMILES string represents a chalcogen atom, including isotopes.
    
    Args:
        smiles (str): SMILES string of a chemical entity

    Returns:
        bool: True if the SMILES represents a chalcogen atom, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Consider atomic numbers for chalcogens: O (8), S (16), Se (34), Te (52), Po (84)
    chalcogen_atomic_numbers = {8, 16, 34, 52, 84}
    
    # Check all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in chalcogen_atomic_numbers:
            # Allow isotopes, which have the same atomic number
            return True, f"The atom is a chalcogen: {atom.GetSymbol()}"
    
    return False, "No chalcogen atom present"