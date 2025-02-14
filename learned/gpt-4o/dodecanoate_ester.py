"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester should have lauric acid (12-carbon saturated chain) as the acid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define lauric acid ester SMARTS pattern
    lauric_acid_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O")
    if lauric_acid_ester_pattern is None:
        return None, None  # Return none if the pattern can't be compiled
    
    # Check for lauric acid ester substructure
    if mol.HasSubstructMatch(lauric_acid_ester_pattern):
        return True, "Contains lauric acid ester component"
    
    return False, "Does not contain lauric acid ester component"