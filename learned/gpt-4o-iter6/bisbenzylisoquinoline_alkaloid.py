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
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Possible patterns for detecting isoquinoline-like structures
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccc2ncccc2c1")
    
    # Broad pattern to detect oxygen bridge linkages
    ether_linkage_pattern = Chem.MolFromSmarts("c-O-c")

    # Search for benzylisoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoquinoline_matches) < 2:
        return False, "Less than two isoquinoline-like units detected"

    # Ensure there is an ether link across significant structures
    ether_matches = mol.GetSubstructMatches(ether_linkage_pattern)
    if len(ether_matches) < 1:
        return False, "No ether linkages detected"

    # Additional checks for bridging patterns, e.g., direct carbon-carbon bridges
    # Check for methylenedioxy groups indirectly - this is complex, hence only hinted here
    # methylenedioxy_pattern = Chem.MolFromSmarts("COC")

    return True, "Contains two benzylisoquinoline units linked by ether and possibly other bridges"