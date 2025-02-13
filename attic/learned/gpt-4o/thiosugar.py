"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has a carbohydrate backbone with one or more oxygens or hydroxy
    groups replaced by sulfur or -SR groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broader patterns to match common sugars and linked hydroxyls (e.g., pyranoses, furanoses)
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O"),  # pyranose
        Chem.MolFromSmarts("C1OC(O)C(O)C1O"),      # furanose
        Chem.MolFromSmarts("O=C(O)C(O)C(O)C"),     # open-chain
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No recognizable carbohydrate motifs found"
    
    # Enhance sulfur substitution detection
    sulfur_substitution_patterns = [
        Chem.MolFromSmarts("[SX2]"),            # Sulfide
        Chem.MolFromSmarts("[SX3]=O"),          # Sulfoxide
        Chem.MolFromSmarts("[SX4](=O)(=O)"),    # Sulfone
        Chem.MolFromSmarts("[SX2]-[C]"),        # Sulfur linked to carbon indicating an -SR group
        Chem.MolFromSmarts("[SH2]"),            # Thiol
        Chem.MolFromSmarts("[SX2H]")            # Thiol form (rare in sugars, but cover it)
    ]
    
    for pattern in sulfur_substitution_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Sulfur substitution found in potential carbohydrate structure"
    
    return False, "No sulfur substitution found in expected carbohydrate structure locations"