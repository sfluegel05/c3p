"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid must contain both a sphingoid/ceramide backbone and a carbohydrate residue
    attached via a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a carbohydrate moiety with a glycosidic linkage pattern
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Look for sphingoid base pattern
    sphingoid_pattern = Chem.MolFromSmarts("O[C@H](CO)CN(C(=O)C)CCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid/ceramide backbone found"

    # Both the patterns should be present for the classification as glycosphingolipid
    return True, "Contains sphingoid/ceramide backbone with glycosidic linkage to carbohydrate"