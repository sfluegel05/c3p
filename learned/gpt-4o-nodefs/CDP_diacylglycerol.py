"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a glycerol backbone, two ester-linked fatty acid chains, and a cytidine diphosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone - C(CO)O as a simple motif
    glycerol_backbone_pattern = Chem.MolFromSmarts("C(CO)O")
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"

    # Look for at least two ester groups (C(=O)O)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for cytidine moiety - broader pattern n1c(C)[nH]c1=O
    cytidine_pattern = Chem.MolFromSmarts("n1cc[nH]c1=O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"
    
    # Look for a diphosphate group
    diphosphate_pattern = Chem.MolFromSmarts("P(O)P(O)")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate group found"

    return True, "Contains glycerol backbone with two ester-linked fatty acid chains and cytidine diphosphate group"