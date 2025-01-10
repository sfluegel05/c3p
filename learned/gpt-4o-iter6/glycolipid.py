"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as containing a lipid portion connected via glycosidic linkage to a sugar moiety.

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
    
    # Enhance lipid backbone patterns
    lipid_backbone_patterns = [
        Chem.MolFromSmarts("C(COC(=O)[CX4])"),  # General glycerolipid ester linkage
        Chem.MolFromSmarts("N[C@@H](C)[CX4]"),  # Amide linkage typical in sphingolipids
        Chem.MolFromSmarts("[CH2]CCCCCCCCCCCCCCC"),  # Long alkane chain (common in lipids)
    ]
    
    if not any(mol.HasSubstructMatch(pat) for pat in lipid_backbone_patterns):
        return False, "No recognizable lipid backbone structure found"

    # Improve glycosidic linkage pattern
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[OX2]C([OX2])[OX2]C")  # Typical glycosidic bond pattern
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage detected"

    # Comprehensive sugar moiety patterns
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC([OX2])[CX4]([OX2])[C@H]([CX4]1)"),  # common pyranose forms
        Chem.MolFromSmarts("C1OC(O)[CH](O)[C@@H]1[C@H]"),  # beta or alpha anomeric
    ]

    if not any(mol.HasSubstructMatch(pat) for pat in sugar_patterns):
        return False, "No recognizable sugar moiety found"
    
    return True, "The structure matches the criteria for a glycolipid, containing a lipid backbone linked via glycosidic linkage to a sugar moiety"