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
    
    # Look for glycerol backbone pattern (C-C-C where two carbons are bound to hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[O][C][C][C](O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for esterified acyl group pattern (O=C-[C] linked to O of glycerol)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C][C][C](O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Count the esterified acyl groups
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} esterified acyl groups, requires exactly 1"

    return True, "Contains a glycerol backbone with one esterified acyl group"