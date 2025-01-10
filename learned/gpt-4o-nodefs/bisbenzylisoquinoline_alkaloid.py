"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    Bisbenzylisoquinoline alkaloids are characterized by two benzylisoquinoline units linked by 
    ether bonds, featuring complex ring systems and methoxy groups.

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
    
    # Define substructure patterns for bisbenzylisoquinoline alkaloids
    ether_linkage_pattern = Chem.MolFromSmarts("CCOc1ccccc1")  # Simple ether linkage and aromatic ring
    methoxy_pattern = Chem.MolFromSmarts("CO")  # Methoxy group
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccc2ncccc2c1")  # Isoquinoline
    
    # Checking for ether linkages
    if not mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "No appropriate ether linkages found"
    
    # Checking for methoxy groups
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 2:
        return False, f"Found {len(methoxy_matches)} methoxy groups, need at least 2"
    
    # Checking for isoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoquinoline_matches) < 2:
        return False, f"Found {len(isoquinoline_matches)} isoquinoline units, need at least 2"

    return True, "Contains characteristic features of bisbenzylisoquinoline alkaloids"