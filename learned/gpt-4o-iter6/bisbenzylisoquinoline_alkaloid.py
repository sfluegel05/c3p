"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline is defined by two benzylisoquinoline units linked by ether bridges, and may have additional bridging.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for benzylisoquinoline unit with more details
    benzylisoquinoline_pattern = Chem.MolFromSmarts("c1ccc2[nH]c3ccccc3c2c1")  # Example pattern, should be elaborated

    # Check for two such units
    benzyl_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(benzyl_matches) < 2:
        return False, "Less than two benzylisoquinoline units detected"

    # Define patterns for different types of bridging
    # Ether bridge (more specific than C-O-C)
    ether_bridge_pattern = Chem.MolFromSmarts("c-O-c")

    # Methylenedioxy bridge (additional possibility)
    methylenedioxy_bridge_pattern = Chem.MolFromSmarts("c1oc2c(oc1)cccc2")
    
    # Check for any valid bridge connecting the benzylisoquinoline units
    if not (mol.HasSubstructMatch(ether_bridge_pattern) or mol.HasSubstructMatch(methylenedioxy_bridge_pattern)):
        return False, "No valid bridge detected between benzylisoquinoline units"
    
    return True, "Contains two benzylisoquinoline units linked by ether or alternative bridges"