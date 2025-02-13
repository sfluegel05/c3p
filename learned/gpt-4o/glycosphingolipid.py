"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid must contain both a sphingoid/ceramide backbone and a carbohydrate residue
    attached via a glycosidic linkage to the O-1 of the sphingoid.

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

    # Look for pattern indicating glycosidic linkage (O-glycosidic bound sugar)
    glycosidic_pattern = Chem.MolFromSmarts("CO[C@H]1O[C@@H]([C@@H](O)[C@@H](O)[C@H]1O)[C@@H]([C@H](O)CO)O")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Look for long-chain base sphingoid/ceramide backbone pattern
    sphingoid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCOC(=O)N[C@@H](CO)C")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid/ceramide backbone found"

    # Both substructures must exist in a glycosphingolipid
    return True, "Contains sphingoid/ceramide backbone with glycosidic linkage to carbohydrate"