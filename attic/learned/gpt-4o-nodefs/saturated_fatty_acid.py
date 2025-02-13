"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carboxylic acid groups (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    num_carboxylic_acid_groups = len(mol.GetSubstructMatches(carboxylic_acid_pattern))
    if num_carboxylic_acid_groups != 1:
        return False, f"Found {num_carboxylic_acid_groups} carboxylic acid groups; requires exactly 1"
    
    # Check for unsaturation: presence of double bonds
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        return False, "Contains unsaturation (double bonds found)"
    
    return True, "Molecule is a saturated fatty acid"