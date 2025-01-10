"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is a glycerol backbone with a phosphate group at position 3,
    and a single acyl group attached at either position 1 or position 2.
    
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
    
    # Neutralize charges (if any)
    Chem.SanitizeMol(mol)

    # Define glycerol 3-phosphate backbone pattern (allowing for stereochemistry)
    glycerol3p_pattern = Chem.MolFromSmarts("""
    [O][P](=O)([O])[O][C@@H]([CH2][OH])[CH2][OH]
    """)
    # Alternatively, account for both enantiomers
    glycerol3p_pattern_rev = Chem.MolFromSmarts("""
    [O][P](=O)([O])[O][C@H]([CH2][OH])[CH2][OH]
    """)

    if not (mol.HasSubstructMatch(glycerol3p_pattern) or mol.HasSubstructMatch(glycerol3p_pattern_rev)):
        return False, "No glycerol 3-phosphate backbone found"
    
    # Check for single acyl group attached via ester linkage at position 1 or 2
    # Define ester linkage pattern with variable attachment point
    ester_pattern = Chem.MolFromSmarts("""
    [C](=O)[O][CH][CH2][O][P]
    """)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Define possible positions for acyl attachment (positions 1 or 2)
    acyl_positions_patterns = [
        Chem.MolFromSmarts("[C](=O)[O][CH]([CH2][OH])[CH2][OH]"),  # Acyl at position 1
        Chem.MolFromSmarts("[C](=O)[O][CH2][CH]([OH])[CH2][OH]"),  # Acyl at position 2
    ]

    acyl_found = False
    for pattern in acyl_positions_patterns:
        if mol.HasSubstructMatch(pattern):
            acyl_found = True
            break

    if not acyl_found:
        return False, "No acyl group attached at position 1 or 2 via ester linkage"
    
    # Check for the number of acyl groups (should be exactly one)
    acyl_group_pattern = Chem.MolFromSmarts("[C](=O)[O][CH]")
    acyl_groups = mol.GetSubstructMatches(acyl_group_pattern)
    if len(acyl_groups) != 1:
        return False, f"Found {len(acyl_groups)} acyl groups, expected exactly 1"
    
    # Ensure there are no additional ester or ether linkages
    ester_linkage_pattern = Chem.MolFromSmarts("[C](=O)[O][CH]")
    ester_linkages = mol.GetSubstructMatches(ester_linkage_pattern)
    if len(ester_linkages) > 1:
        return False, "Additional ester linkages found"

    ether_linkage_pattern = Chem.MolFromSmarts("[C][O][CH]")
    ether_linkages = mol.GetSubstructMatches(ether_linkage_pattern)
    if len(ether_linkages) > 0:
        return False, "Ether linkages found"

    # Confirm that the phosphate group is at position 3
    phosphate_pattern = Chem.MolFromSmarts("[C][CH2][O][P](=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not attached at position 3"

    return True, "Contains glycerol 3-phosphate backbone with a single acyl group attached at position 1 or 2"