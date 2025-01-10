"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol moiety pattern (C-C-C where two have hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[CX4](CO)[CX4](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for esterified acyl group pattern (O=C-C linked to O of glycerol)
    ester_pattern = Chem.MolFromSmarts("O=C([#6])O[CX4](CO)[CX4](CO)CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No esterified acyl group linked to glycerol found"
    
    # Count the esterified acyl groups
    ester_groups = mol.GetSubstructMatches(ester_pattern)
    if len(ester_groups) != 1:
        return False, f"Found {len(ester_groups)} esterified acyl groups, requires exactly 1"
    
    # Verify the molecular structure: ensure only one esterified acyl group
    if len(ester_groups) != 1:
        return False, "More than one esterified acyl group"
    
    return True, "Contains a glycerol backbone with one esterified acyl group"