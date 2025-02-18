"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid must contain both a sphingoid/ceramide backbone
    and a carbohydrate residue attached via a glycosidic linkage to O-1 of the sphingoid.
    
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

    # Refined general glycosidic linkage pattern
    glycosidic_pattern = Chem.MolFromSmarts("CO[C;R]")  # Oxygen linking two carbons, C in a ring
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Refined sphingoid base, which includes long chain and amide link
    sphingoid_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO[C])C(O)[C;!R][C;!R]")  
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid/ceramide backbone found"

    # Both the patterns should be present for the classification as glycosphingolipid
    return True, "Contains sphingoid/ceramide backbone with glycosidic linkage to carbohydrate"