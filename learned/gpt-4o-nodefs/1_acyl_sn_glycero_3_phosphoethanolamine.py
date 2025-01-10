"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This is identified by a chiral glycerol backbone, a single acyl group, and phosphoethanolamine.

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
    # Must have a chiral center on the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC(=O)C)[CH2]OP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No chiral glycerol backbone found with the specific arrangement"

    # Acyl group pattern - count occurrence to ensure only one
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Incorrect number of acyl groups found: expected 1, found {len(acyl_matches)}"

    # Phosphoethanolamine pattern
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("OP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phospho_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # If all features are present, classify as 1-acyl-sn-glycero-3-phosphoethanolamine
    return True, "Contains chiral glycerol backbone, single acyl group, and phosphoethanolamine group"