"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This is identified by a chiral glycerol backbone, acyl group, and phosphoethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Chiral glycerol backbone pattern
    # This pattern specifies the chiral sn-glycerol with one specific 2R-hydroxy position.
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COC(=O)[C])[CH2]OP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No chiral glycerol backbone found"

    # Acyl group typically recognized as an ester linkage on the glycerol
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached through ester linkage found"

    # Phosphoethanolamine pattern accounting for common phospholipid linkage and possible protonation states
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("OP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phospho_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # If all features are present, classify as 1-acyl-sn-glycero-3-phosphoethanolamine
    return True, "Contains chiral glycerol backbone, acyl group, and phosphoethanolamine group"