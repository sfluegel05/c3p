"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline is defined by two benzylisoquinoline units linked by ether bridges, with possible additional bridging.

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

    # Define SMARTS for bisbenzylisoquinoline pattern
    # This pattern looks for isoquinoline units with attached ether linkages
    isoquinoline_unit_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cc[nH]c2")
    ether_bridge_pattern = Chem.MolFromSmarts("-O-")

    # Check for two isoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_unit_pattern)
    if len(isoquinoline_matches) < 2:
        return False, "Less than two isoquinoline units detected"

    # Check for ether bridge between isoquinoline units
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if not ether_matches:
        return False, "No valid ether bridges detected between units"

    # Additional bridging check (optional, based on complex structures seen)
    methylenedioxy_bridge_pattern = Chem.MolFromSmarts("O1COC2=CC(=O)NC2")
    if mol.HasSubstructMatch(methylenedioxy_bridge_pattern):
        return True, "Contains two isoquinoline units linked by ether bridges and additional methylenedioxy bridging"

    return True, "Contains two isoquinoline units linked by ether bridges"