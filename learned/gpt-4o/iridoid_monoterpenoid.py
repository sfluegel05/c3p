"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid contains a cyclopentane ring fused to a six-membered oxygen heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for iridoid and secoiridoid structures
    iridoid_pattern = Chem.MolFromSmarts("C1C=CC2OC2C1")  # Generic pattern for cyclopentane and six-membered ether ring
    secoiridoid_pattern = Chem.MolFromSmarts("C1CC(C=O)C2OC2C1")  # Pattern for secoiridoids with typical ring openings

    # Check for iridoid and secoiridoid structures
    is_iridoid = mol.HasSubstructMatch(iridoid_pattern)
    is_secoiridoid = mol.HasSubstructMatch(secoiridoid_pattern)
    
    if is_iridoid:
        return True, "Contains an iridoid structure with a fused cyclopentane and six-membered oxygen heterocycle."
    
    if is_secoiridoid:
        return True, "Contains a secoiridoid structure with a characteristic ring opening."
    
    return False, "Does not match typical iridoid monoterpenoid structures"