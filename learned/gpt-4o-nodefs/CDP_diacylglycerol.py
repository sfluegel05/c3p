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
        bool: True if the molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined glycerol backbone pattern (OCC(COC)OC)
    glycerol_backbone_pattern = Chem.MolFromSmarts("OCC(COC)OC")
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"

    # Look for at least two ester groups (C(=O)OC)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Improved cytidine moiety pattern
    cytidine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c([nH]1)=O")
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"
    
    # Look for a diphosphate group - P(~[OH])([OH])(~[O])~O
    diphosphate_pattern = Chem.MolFromSmarts("P([OH])([OH])([O])~[O]")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate group found"

    return True, "Contains glycerol backbone with two ester-linked fatty acid chains and cytidine diphosphate group"