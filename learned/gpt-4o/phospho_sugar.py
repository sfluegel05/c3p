"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy 
    group esterified with phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General sugar backbone pattern
    # Simplified to focus on serial CO pattern indicative of pendulous saccharide chains
    sugar_pattern = Chem.MolFromSmarts("[C&H][O&H][C&H][O&H]")  # More flexible for sugar backbone
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like structure found"
    
    # General phosphate ester linkage
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")  # Phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester linkage found"
    
    return True, "Contains a saccharide with a phosphate ester linkage"