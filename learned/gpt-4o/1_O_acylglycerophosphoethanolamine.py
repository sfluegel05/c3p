"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for 1-O-acylglycerophosphoethanolamine components
    
    # 1-O-acylglycerol pattern: more flexible stereochemistry and backbone chains
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COC(=O)[C])[CH2]OC")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1-O-acyl-glycerol structure found"

    # Phosphoethanolamine group pattern allowing variations in charge representation
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not present"
    
    # Check if components link correctly to match a 1-O-acylglycerophosphoethanolamine
    if not (mol.HasSubstructMatch(glycerol_pattern) and mol.HasSubstructMatch(phosphoethanolamine_pattern)):
        return False, "Composite structure not consistent with 1-O-acylglycerophosphoethanolamine"
    
    return True, "Structure matches 1-O-acylglycerophosphoethanolamine"