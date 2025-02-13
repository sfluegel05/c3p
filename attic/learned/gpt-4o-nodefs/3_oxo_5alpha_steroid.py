"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for 3-oxo-5alpha-steroid.
    # This pattern includes a steroid backbone with a 3-oxo group.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(=O)CC4C(C3=CC2C1)CCC4")    
    if mol.HasSubstructMatch(steroid_pattern):
        # Check if the molecule has the steroid backbone pattern
        if mol.HasSubstructMatch(steroid_pattern):
            return True, "Matches the core 3-oxo-5alpha-steroid pattern"
        else:
            return False, "Does not match the steroid backbone pattern"
    else:
        return False, "Does not contain 3-oxo group at required position"

    return None, None