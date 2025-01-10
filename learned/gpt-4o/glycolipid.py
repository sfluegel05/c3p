"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a type of lipid that includes a carbohydrate group and usually consists 
    of a 1,2-di-O-acylglycerol unit linked at the oxygen 3 by a glycosidic linkage to a carbohydrate.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol with acyl chains pattern
    glycerol_acyl_pattern = Chem.MolFromSmarts("C(CO*)O[*]C(C=O)O")  # Simplified SMARTS for diacylglycerol
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "No diacylglycerol structure found"
        
    # Look for carbohydrate (sugar) moiety, common patterns for monosaccharides
    sugar_patterns = ["O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O",  # Glucose-like pattern
                      "O[C@H]1[C@H](O)[C@H](O)[C@@H](CO)O[C@@H]1"]  # Galactose-like pattern
    sugar_found = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No carbohydrate moiety found"

    return True, "Contains a 1,2-di-O-acylglycerol unit linked to a carbohydrate"