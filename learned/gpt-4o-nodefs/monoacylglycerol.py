"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Look for glycerol backbone pattern (Three carbon atoms each with at least one oxygen)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)OCC(O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for one ester group pattern (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # All checks passed, the molecule is a monoacylglycerol
    return True, "Contains glycerol backbone with one esterified fatty acid chain"