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
    A bisbenzylisoquinoline is defined by two benzylisoquinoline units linked by ether bridges.

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
    
    # Define SMARTS for benzylisoquinoline unit (simplified)
    benzylisoquinoline_pattern = Chem.MolFromSmarts("c1ccccc1Cc1cc2c(ccc3CCN(CCc4ccccc4)c3)c2cc1")

    # Check for two such units connected by oxygen linkage
    ether_bridge_pattern = Chem.MolFromSmarts("c-O-c")
    
    # Check for two benzylisoquinoline units
    if len(mol.GetSubstructMatches(benzylisoquinoline_pattern)) < 2:
        return False, "Less than two benzylisoquinoline units detected"
    
    # Check for ether bridges
    if not mol.HasSubstructMatch(ether_bridge_pattern):
        return False, "No ether bridge detected between units"

    return True, "Contains two benzylisoquinoline units linked by ether bridges"