"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol contains a glycerol backbone where one position is an esterified acyl group
    and the remaining two are either H or alkyl groups.

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
    
    # Look for glycerol backbone pattern (C-C-C with two hydroxyls)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for esterified acyl group (OC(=O) chain)
    ester_pattern = Chem.MolFromSmarts("OC(=O)[CX4,CX3]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check that exactly one esterified acyl group is present
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} esterified acyl groups, needs exactly 1"

    return True, "Contains a glycerol backbone with one esterified acyl group"