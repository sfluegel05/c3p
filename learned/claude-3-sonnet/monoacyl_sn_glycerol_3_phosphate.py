"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:57980 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
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

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for ester group (acyl chain)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for glycerol backbone with correct substitution pattern
    # [OH]-C-C(OH)-C-O-P pattern or [OH]-C-C(O-acyl)-C-O-P pattern
    glycerol_phosphate_pattern1 = Chem.MolFromSmarts("[OX2H]C[CH]([OX2H,OX2C(=O)])[CH2]OP(=O)([OX2H,OX1-])[OX2H,OX1-]")
    glycerol_phosphate_pattern2 = Chem.MolFromSmarts("[OX2H,OX2C(=O)]C[CH]([OX2H])[CH2]OP(=O)([OX2H,OX1-])[OX2H,OX1-]")
    
    if not (mol.HasSubstructMatch(glycerol_phosphate_pattern1) or mol.HasSubstructMatch(glycerol_phosphate_pattern2)):
        return False, "No glycerol backbone with correct substitution pattern found"

    # Count carbons in acyl chain (should be at least 4)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2][CH]([OX2H,OX2C(=O)])[CH2]OP")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No proper acyl chain attachment found"

    # Verify only one acyl chain by counting ester carbonyls
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if len(mol.GetSubstructMatches(carbonyl_pattern)) != 1:
        return False, "Incorrect number of acyl groups"

    # Count free hydroxyls on glycerol backbone
    free_oh_pattern = Chem.MolFromSmarts("[CH2,CH]([CH2,CH][OX2H,OX2P,OX2C(=O)])[OX2H]")
    if not mol.HasSubstructMatch(free_oh_pattern):
        return False, "Missing free hydroxyl group on glycerol backbone"

    return True, "Contains glycerol backbone with one acyl chain and phosphate group at position 3"