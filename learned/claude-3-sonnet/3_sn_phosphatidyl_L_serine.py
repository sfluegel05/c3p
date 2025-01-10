"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:57262 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has:
    - A glycerol backbone with R stereochemistry at sn-2
    - Two acyl groups at sn-1 and sn-2 positions
    - A phosphoserine group at sn-3 with L-serine stereochemistry
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check: Complete core structure with correct stereochemistry
    # This SMARTS pattern defines:
    # - Glycerol backbone with R stereochemistry at sn-2 [C@H]
    # - Two ester groups at sn-1 and sn-2
    # - Phosphate at sn-3
    # - L-serine with S stereochemistry [C@H]
    core_pattern = Chem.MolFromSmarts(
        "[CH2X4][C@H]([CH2X4][OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O)[OX2]"
    )
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core structure or stereochemistry incorrect"

    # Check for exactly two ester groups
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2,CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for exactly one phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=[O])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, "Must have exactly one phosphate group"

    # Verify L-serine moiety
    serine_pattern = Chem.MolFromSmarts("[NH2][C@H](C(=O)[OH])CO")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No L-serine moiety found or incorrect stereochemistry"

    # Verify acyl groups have proper fatty acid chains
    # Looking for chains of at least 4 carbons
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing proper fatty acid chains"

    # Additional checks for atom counts
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if p_count != 1 or n_count != 1:
        return False, "Incorrect number of phosphorus or nitrogen atoms"

    # More specific connectivity check for phosphoserine at sn-3
    sn3_pattern = Chem.MolFromSmarts(
        "[CH2X4][C@H]([OX2]C(=O)*)([CH2X4][OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O)"
    )
    if not mol.HasSubstructMatch(sn3_pattern):
        return False, "Incorrect connectivity at sn-3 position"

    return True, "Contains glycerol backbone with two acyl groups and phosphoserine moiety in correct configuration"