"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a steroid backbone pattern (tetracyclic core)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Look for glycosidic linkage patterns (C-O-C linkage indicative of sugars)
    glycosidic_pattern = Chem.MolFromSmarts("[C,O]-O-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic linkage found"
        
    # Check for typical sugar moiety (e.g., glucose, based on OH groups and specific patterns)
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C(CO)O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moiety found"

    return True, "Contains steroid backbone with glycosidic linkage to sugar moieties"