"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined by the presence of a 1,2-di-O-acylglycerol structure 
    linked to a glycosidic saccharide moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Generic glycerol backbone with acyl groups pattern
    glycerol_diacyl_pattern = Chem.MolFromSmarts("C(COC(=O)[CX4][O])OC(=O)[CX4]")
    if not mol.HasSubstructMatch(glycerol_diacyl_pattern):
        return False, "No 1,2-di-O-acylglycerol structure found"

    # Generic glycosidic linkage to a carbohydrate moiety
    glycosidic_linkage_pattern = Chem.MolFromSmarts("OC[C@]([AX4])(O[C@]1([AX4])O[C@H]([AX4])C([AX4])O1)")
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage with carbohydrate moiety detected"

    # Enhanced pattern for recognizing a sugar moiety (broad coverage)
    sugar_pattern = Chem.MolFromSmarts("C1OC(C[OX2H])C[OX2H]C1[OX2H]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognizable sugar moiety found"
    
    return True, "The structure matches the criteria for a glycolipid with 1,2-di-O-acylglycerol structure and glycosidic linkage to sugar moiety"