"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: CHEBI:flavin derivatives (dimethylisoalloxazine derivatives)
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    Flavin must contain the dimethylisoalloxazine core with a substituent at position 10.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern using SMARTS
    # Pattern matches dimethylisoalloxazine core with substituent at N10
    core_smarts = '[CH3]c1c([CH3])cc2nc3c(=O)[nH]c(=O)nc3n([!H0])c12'
    core_pattern = Chem.MolFromSmarts(core_smarts)
    
    # Check for core structure match
    if mol.HasSubstructMatch(core_pattern):
        return True, "Contains dimethylisoalloxazine core with substituent at position 10"
    
    # Additional check for possible variations (e.g., different oxidation states)
    # Match core without explicit methyl hydrogens but with methyl substitution
    alt_smarts = 'c1(C)c(C)cc2nc3c(=O)[nH]c(=O)nc3n([!H0])c12'
    alt_pattern = Chem.MolFromSmarts(alt_smarts)
    
    if mol.HasSubstructMatch(alt_pattern):
        return True, "Contains dimethylisoalloxazine core with substituent at position 10"

    return False, "Missing dimethylisoalloxazine core or position 10 substituent"