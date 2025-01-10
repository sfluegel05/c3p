"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A molecule is classified as such if it contains a glycerol backbone with correct stereochemistry,
    acyl groups at the sn-1 and sn-2 positions, and a phosphoserine group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with sn-3 stereochemistry
    # Modify SMARTS to reflect accurate stereochemistry
    glycerol_sn3_pattern = Chem.MolFromSmarts("O[C@H](COP(O)(=O)O)[C@@H](O)")
    if not mol.HasSubstructMatch(glycerol_sn3_pattern):
        return False, "No proper glycerol sn-3 stereochemistry found"

    # Ensure two acyl groups on sn-1 and sn-2
    # Recognize ester linkages, confirm correct stereochemistry
    acyl_pattern_sn1_sn2 = Chem.MolFromSmarts("C(=O)O[C@H]")
    acyl_matches_sn1_sn2 = mol.GetSubstructMatches(acyl_pattern_sn1_sn2)
    if len(acyl_matches_sn1_sn2) < 2:
        return False, f"Found {len(acyl_matches_sn1_sn2)} acyl groups, need exactly 2"

    # Check for the phosphoserine group at sn-3
    phosphoserine_pattern = Chem.MolFromSmarts("COP(O)(=O)O[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # If all patterns match, classify as 3-sn-phosphatidyl-L-serine
    return True, "Contains the structures indicative of a 3-sn-phosphatidyl-L-serine"