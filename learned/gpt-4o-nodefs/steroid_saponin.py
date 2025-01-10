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

    # Flexible steroid backbone pattern including stereochemistry
    steroid_pattern = Chem.MolFromSmarts("C1CC2=CC3C=C[C@H](CC3)C2(C)CC4C1(C)CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Enhanced glycosidic linkage pattern allowing for variable sugar attachments
    glycosidic_pattern = Chem.MolFromSmarts("[C,O]-[O]-[C,O]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic linkage found"
        
    # Broader sugar moiety recognition, allowing for common sugar patterns
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)C")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moiety found"

    return True, "Contains steroid backbone with glycosidic linkage to sugar moieties"