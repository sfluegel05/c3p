"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline is defined by two benzylisoquinoline units linked by ether bridges,
    with additional possible bridging such as carbon-carbon bridges or methylenedioxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted pattern for detecting isoquinoline-like structures
    # Allow benzyl isoquinoline variations by detecting segments typical for bisbenzylisoquinoline
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccc2ncccc2c1")

    # Improved pattern to detect ether linkages, including aromatic oxygens
    ether_linkage_pattern = Chem.MolFromSmarts("cOc")

    # Search for benzylisoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoquinoline_matches) < 2:
        return False, "Less than two benzylisoquinoline-like units detected"

    # Ensure there is at least one ether linkage
    ether_matches = mol.GetSubstructMatches(ether_linkage_pattern)
    if len(ether_matches) < 1:
        return False, "No ether linkages detected"

    # Check for additional bridging patterns, such as methylenedioxy groups
    methylenedioxy_pattern = Chem.MolFromSmarts("OCcO")
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    
    # Consider either ether or methylenedioxy bridging for at least one linking
    if (len(ether_matches) + len(methylenedioxy_matches)) < 1:
        return False, "Missing bridging between two units"

    return True, "Contains two benzylisoquinoline units linked by ether bridges and possible additional bridging"