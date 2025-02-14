"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol in which the acyl substituent is located at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone"

    # Check for single ester linkage attached to the second carbon
    ester_pattern = Chem.MolFromSmarts("OCC(OC(=O)[#6])CO") # Acyl group off second carbon
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Ester group not on the second carbon"

    # Verify that there is only one ester group
    general_ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(general_ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Ensure no additional esterification occurs at positions other than the 2-position
    non_2_position_ester_pattern = Chem.MolFromSmarts("[OX2][#6][OX2H]")
    other_ester_matches = mol.GetSubstructMatches(non_2_position_ester_pattern)
    if len(other_ester_matches) > 1:
        return False, f"Ester group is not exclusively on the 2-position"

    return True, "Contains a glycerol backbone with an acyl substituent at position 2"