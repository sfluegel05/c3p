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
    
    # 1-O-acylglycerol pattern: stereochemistry for the glycerol backbone
    # [C@@H] for stereochemistry, the central carbon with oxygens on each side
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC(=O)[C])[CH2]OC")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1-O-acyl-glycerol structure found"

    # Phosphoethanolamine group pattern
    # Ensure that phosphate links through the correct sites with ethanolamine 
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCC[NH2+]")  
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not present"
    
    # Assess complete molecule structure for exact match to 1-O-acylglycerophosphoethanolamine description
    if not (mol.HasSubstructMatch(glycerol_pattern) and mol.HasSubstructMatch(phosphoethanolamine_pattern)):
        return False, "Composite structure not consistent with 1-O-acylglycerophosphoethanolamine"
    
    return True, "Structure matches 1-O-acylglycerophosphoethanolamine"