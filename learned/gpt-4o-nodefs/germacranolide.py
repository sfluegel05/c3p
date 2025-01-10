"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are a subclass of sesquiterpene lactones typically having
    a cyclodecadiene core with a lactone group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broadening lactone pattern check
    # Check for 5- and 6-membered lactone rings with possible saturation/unsaturation and common variations
    lactone_patterns = [
        Chem.MolFromSmarts("[C;R][O][C;R](=O)"),  # General lactone pattern
        Chem.MolFromSmarts("O[C;R]=C"),           # Alternative oxygen and carbon positions
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "Lactone group not found"
    
    # Check for decalin (fused cyclohexane) systems as germacranolides typically have complex ring structures
    decalin_pattern = Chem.MolFromSmarts("C1CC2CCCC(C1)C2")  # Generalized large ring system, often found in germacranolides
    
    if not mol.HasSubstructMatch(decalin_pattern):
        return False, "No decalin or cyclodecadiene structure found"
    
    # Check for multiple double bonds in ring systems
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 3:  # Increasing threshold for more specificity
        return False, "Insufficient number of double bonds for germacranolides"
    
    return True, "Contains a characteristic lactone group and complex ring system of germacranolides"