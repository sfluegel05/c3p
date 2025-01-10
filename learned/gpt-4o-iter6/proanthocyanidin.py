"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer consisting of two or more hydroxyflavan units, often linked by C4-C8 or C4-C6 bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for hydroxy groups in phenolic rings as part of flavan units
    phenol_pattern = Chem.MolFromSmarts("c1c(O)ccc(O)c1")  # Phenolic pattern
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) < 4:  # Ensure at least indicative of two flavan units
        return False, "Insufficient phenolic groups to suggest multiple hydroxyflavan units"
    
    # Patterns for interflavan C4-C8 or C4-C6 linkages
    linkage_pattern_c8 = Chem.MolFromSmarts("[C@H]1c(O)cc(O)c1-c2c(O)ccc(O)c2")  # One approach to C4-C8 linkage with stereochemistry
    linkage_pattern_c6 = Chem.MolFromSmarts("[C@H]1c(O)cc(O)c1-c2cc(O)ccc2")  # One approach to C4-C6 linkage
    
    # Check for presence of any valid linkage
    linkage_matches_c8 = mol.GetSubstructMatches(linkage_pattern_c8)
    linkage_matches_c6 = mol.GetSubstructMatches(linkage_pattern_c6)
    
    if len(linkage_matches_c8) < 1 and len(linkage_matches_c6) < 1:
        return False, "No interflavan linkage found"
    
    return True, "Contains multiple hydroxyflavan units with at least one interflavan linkage"