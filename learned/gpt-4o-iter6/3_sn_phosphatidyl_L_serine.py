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

    # Glycerol backbone recognition with sn-3 stereochemistry
    # [C@H](CO) pattern for stereospecificity, linking all critical groups
    glycerol_sn_pattern = Chem.MolFromSmarts("O[C@@H](COP(=O)(O)OC[C@H](N)C(=O)O)[C@@H](OC(=O))")
    if not mol.HasSubstructMatch(glycerol_sn_pattern):
        return False, "No proper glycerol sn-3 stereochemistry found"

    # Acyl groups on sn-1 and sn-2
    # Recognize ester linkages - virtually "C(=O)O[C@H]"
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 2:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 2"

    # Phosphoserine moiety recognition
    # Looking specifically for phosphate group linked to serine
    phosphoserine_pattern = Chem.MolFromSmarts("COP(=O)(OC[C@H](N)C(=O)O)O")
    phosphoserine_match = mol.HasSubstructMatch(phosphoserine_pattern)
    if not phosphoserine_match:
        return False, "No phosphoserine group found"

    # If all patterns match, classify as 3-sn-phosphatidyl-L-serine
    return True, "Contains the structures indicative of a 3-sn-phosphatidyl-L-serine"