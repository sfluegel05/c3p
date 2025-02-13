"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.

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

    # Look for glycerol backbone pattern (C-C-C backbone)
    glycerol_pattern = Chem.MolFromSmarts("CC(C)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for choline part of phosphocholine group ([N+](C)(C)C)
    choline_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(choline_pattern):
        return False, "No choline group found"

    # Look for phosphate group attached to choline (P(=O)([O-])OC)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])OC")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for two ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    return True, "Contains glycerol backbone with phosphocholine group and 2 esterified fatty acid chains"