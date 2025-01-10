"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glycerol backbone with one esterified fatty acid chain.

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
    
    # Adjust pattern to recognize glycerol backbone more flexibly
    # Use a general pattern for a three-carbon, two-oxygen motif in glycerol
    glycerol_pattern = Chem.MolFromSmarts("[CX4;H2][CX4;H1]([OX2H])[CX4;H2,OX2H]")  # Covers sn substitutions
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for one ester group pattern (C(=O)O), associated with glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # All checks passed, the molecule is a monoacylglycerol
    return True, "Contains glycerol backbone with one esterified fatty acid chain"