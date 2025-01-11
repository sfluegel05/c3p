"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride has a glycerol backbone with two ester-linked fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with at least two oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("C(=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2 or len(ester_matches) > 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    return True, "Contains glycerol backbone with two ester-linked fatty acid chains"