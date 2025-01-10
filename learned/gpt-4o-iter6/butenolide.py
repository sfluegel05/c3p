"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide consists of a 2-furanone skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a 2-furanone core allowing for substitutions
    # The pattern captures the core structure of a gamma-lactone 2-furanone
    furanone_pattern = Chem.MolFromSmarts("C1=CO[C@@H](=O)C1")
    
    # Match furanone skeleton with substitutions
    if mol.HasSubstructMatch(furanone_pattern):
        return True, "Contains 2-furanone skeleton with allowed substitutions"
    else:
        return False, "Does not contain 2-furanone skeleton or allowed substitutions"