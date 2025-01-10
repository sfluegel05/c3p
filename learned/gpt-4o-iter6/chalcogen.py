"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    Chalcogens are elements belonging to group 16 of the periodic table, including O, S, Se, Te, Po.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chalcogen elements
    chalcogen_elements = {'O', 'S', 'Se', 'Te', 'Po'}

    # Check if any atom in the molecule is a chalcogen
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in chalcogen_elements:
            return True, f"Contains chalcogen element: {symbol}"
    
    return False, "No chalcogen elements found"