"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:79126 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has one or more oxygen/hydroxyl groups replaced by sulfur or -SR groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
        
    # Basic carbohydrate check: look for pyranose/furanose ring with multiple hydroxyls
    pyranose_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)S1")
    furanose_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)S1")
    if not (mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)):
        return False, "No carbohydrate backbone detected"
    
    # Check for sulfur substitutions (S attached to sugar backbone)
    sulfur_pattern = Chem.MolFromSmarts("[C,S]-[SX2]")
    sulfur_matches = mol.GetSubstructMatches(sulfur_pattern)
    if not sulfur_matches:
        return False, "No sulfur substitutions found"
    
    # Check for thiol groups (-SH) or thioethers (-S-R)
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    thioether_pattern = Chem.MolFromSmarts("[SX2]-[#6]")
    if not (mol.HasSubstructMatch(thiol_pattern) or mol.HasSubstructMatch(thioether_pattern)):
        return False, "No thiol/thioether groups present"
    
    return True, "Sulfur substitution detected in carbohydrate structure"