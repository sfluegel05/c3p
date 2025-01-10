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

    # Check for exact 3-sn configuration with correct stereochemistry
    # This pattern specifically checks for:
    # - sn-1 ester linkage
    # - sn-2 position with R stereochemistry
    # - sn-3 phosphate linkage
    # - L-serine with S stereochemistry
    core_pattern = Chem.MolFromSmarts(
        "[CH2X4][OX2]C(=O)*.[C@H]([CH2X4][OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O)[OX2]C(=O)*"
    )
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core structure or stereochemistry incorrect"

    # Verify presence of two ester groups
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2,CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check phosphoserine linkage
    phosphoserine_pattern = Chem.MolFromSmarts(
        "[OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O"
    )
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Incorrect phosphoserine group or stereochemistry"

    # Verify fatty acid chains (allowing for modifications)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2][C@H]")
    if len(mol.GetSubstructMatches(fatty_acid_pattern)) < 2:
        return False, "Missing proper acyl groups at sn-1 and sn-2 positions"

    # Check phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=[O])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, "Must have exactly one phosphate group"

    # Verify complete connectivity
    sn3_pattern = Chem.MolFromSmarts(
        "[CH2X4][C@H]([OX2]C(=O)*)([CH2X4][OX2]P(=[O])([OX2])[OX2]C[C@H](N)C(=O)O)"
    )
    if not mol.HasSubstructMatch(sn3_pattern):
        return False, "Incorrect connectivity at sn-3 position"

    # Check for minimum carbon chain length in acyl groups
    chain_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, "Acyl chains too short"

    # Check phosphorus (but not nitrogen to allow for modifications)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"

    return True, "Contains glycerol backbone with two acyl groups and phosphoserine moiety in correct configuration"