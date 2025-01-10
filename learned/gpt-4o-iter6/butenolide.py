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
        bool: True if the molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for a 2-furanone core with substitution tolerance
    # Represents a 5-membered lactone ring with variability in substitution on the ring atoms
    furanone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")  # More generalized to capture various substitutions
    
    # Match furanone skeleton with possible substitutions
    if mol.HasSubstructMatch(furanone_pattern):
        return True, "Contains 2-furanone skeleton with allowed substitutions"
    else:
        return False, "Does not contain 2-furanone skeleton or allowed substitutions"