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

    # Define SMARTS for benzylisoquinoline unit (this should capture common features of isoquinolines)
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccccc2c1ccn2")  # Captures the isoquinoline core
    benzyl_substituent_pattern = Chem.MolFromSmarts("c1ccccc1")  # Simple benzyl group

    # Check for two such units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    benzyl_matches = mol.GetSubstructMatches(benzyl_substituent_pattern)

    if len(isoquinoline_matches) < 2 or len(benzyl_matches) < 2:
        return False, "Less than two benzylisoquinoline units detected"

    # Define patterns for ether bridging
    ether_bridge_pattern = Chem.MolFromSmarts("c-O-c")

    # Check for any valid bridge connecting the benzylisoquinoline units
    if not mol.HasSubstructMatch(ether_bridge_pattern):
        return False, "No valid ether bridge detected between benzylisoquinoline units"

    return True, "Contains two benzylisoquinoline units linked by ether bridges"