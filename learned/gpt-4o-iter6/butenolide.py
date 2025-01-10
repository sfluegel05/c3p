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

    # Define a general SMARTS pattern for a 2-furanone core with substitution tolerance
    # Pattern updated to capture a 5-membered lactone with a double-bond and keto group
    # Allow potential substitutions around the ring structure
    furanone_pattern = Chem.MolFromSmarts("O=C1C=COC1")

    # Match furanone skeleton with allowed substitutions
    if mol.HasSubstructMatch(furanone_pattern):
        return True, "Contains 2-furanone skeleton with allowed substitutions"
    else:
        return False, "Does not contain 2-furanone skeleton or allowed substitutions"