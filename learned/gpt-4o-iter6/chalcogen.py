"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen element based on its SMILES string.
    Chalcogens are elements belonging to group 16 of the periodic table, including O, S, Se, Te, Po.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a single chalcogen atom or isotope, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "SMILES does not represent a single atom"

    # Define chalcogen elements
    chalcogen_elements = {'O', 'S', 'Se', 'Te', 'Po'}

    # Check if the single atom in the molecule is a chalcogen
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()
    if symbol in chalcogen_elements:
        return True, f"Contains chalcogen element: {symbol}"
    
    return False, "No chalcogen elements found"