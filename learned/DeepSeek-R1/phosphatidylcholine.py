"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: CHEBI:16347 phosphatidylcholine
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine has a glycerol backbone with two acyl chains (ester groups)
    and a phosphocholine group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern using SMARTS
    # Pattern matches glycerol backbone with:
    # - Two ester groups at positions 1 and 2
    # - Phosphocholine group at position 3
    core_pattern = Chem.MolFromSmarts(
        "[CH2](-O-C(=O)-*)[CH](-O-C(=O)-*)[CH2]-O-P(=O)([O-])-O-CC[N+](C)(C)C"
    )

    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core phosphatidylcholine structure"

    # Check for exactly two ester groups (acyl chains)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    return True, "Contains glycerol backbone with two acyl chains and phosphocholine group"