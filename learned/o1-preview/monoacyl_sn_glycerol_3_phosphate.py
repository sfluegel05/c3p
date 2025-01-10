"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is a glycerol backbone with a phosphate group at position 3,
    and a single acyl group attached via ester linkage at either position 1 or position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is sanitized
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        return False, "Molecule could not be sanitized"

    # Define glycerol 3-phosphate backbone pattern
    # [C@@H](O)[CH2]O[P](=O)(O)O - represents the sn-glycerol-3-phosphate backbone
    glycerol3p_smarts = "[C@@H](O)[CH2]O[P](=O)(O)O"
    glycerol3p_pattern = Chem.MolFromSmarts(glycerol3p_smarts)
    if glycerol3p_pattern is None:
        return False, "Invalid glycerol 3-phosphate SMARTS pattern"

    # Check for glycerol 3-phosphate backbone
    if not mol.HasSubstructMatch(glycerol3p_pattern):
        return False, "No glycerol 3-phosphate backbone found"

    # Define acyl group attached via ester linkage at position 1 or 2
    # Ester linkage pattern: [C](=O)O[C@@H] or [C](=O)O[CH2]
    ester_pattern1 = Chem.MolFromSmarts("[C](=O)O[C@@H]")
    ester_pattern2 = Chem.MolFromSmarts("[C](=O)O[CH2]")
    if ester_pattern1 is None or ester_pattern2 is None:
        return False, "Invalid ester linkage SMARTS pattern"

    # Check for acyl group attached at position 1 or 2
    acyl_matches1 = mol.GetSubstructMatches(ester_pattern1)
    acyl_matches2 = mol.GetSubstructMatches(ester_pattern2)
    total_acyl_groups = len(acyl_matches1) + len(acyl_matches2)

    if total_acyl_groups != 1:
        return False, f"Expected exactly 1 acyl group, found {total_acyl_groups}"

    # Ensure there are no additional acyl groups
    # Define general ester linkage pattern
    ester_general_pattern = Chem.MolFromSmarts("[C](=O)O[C]")
    ester_general_matches = mol.GetSubstructMatches(ester_general_pattern)
    if len(ester_general_matches) > 1:
        return False, "Additional ester linkages found"

    # Ensure no ether linkages are present
    ether_pattern = Chem.MolFromSmarts("[C]O[C]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) > 0:
        return False, "Ether linkages found"

    return True, "Molecule is a monoacyl-sn-glycerol 3-phosphate"