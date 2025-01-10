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

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone with specific stereochemistry at sn-2 position
    # [C@H] indicates R stereochemistry at carbon 2 of glycerol (sn notation)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H]([CH2X4])[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry found"

    # Check for L-serine moiety with correct stereochemistry
    # [C@H] indicates S stereochemistry for L-serine
    serine_pattern = Chem.MolFromSmarts("[NH2][C@H](C(=O)[OH])CO")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No L-serine moiety found with correct stereochemistry"

    # Check for two ester groups at sn-1 and sn-2 positions
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2,CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Verify the complete 3-sn-phosphatidylserine structure
    # This pattern checks the entire core structure with correct stereochemistry
    complete_pattern = Chem.MolFromSmarts(
        "[CH2X4][C@H]([CH2X4][OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O)[OX2]"
    )
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "Incorrect connectivity or stereochemistry in core structure"

    # Check for phosphorus count
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"

    # Verify fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Fatty acid chains too short or missing"

    # Check that both ester groups are connected to fatty acid chains
    acyl_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, "Missing acyl groups with proper fatty acid chains"

    return True, "Contains glycerol backbone with two acyl groups and phosphoserine moiety in correct configuration"