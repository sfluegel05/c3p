"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol consists of a glycerol backbone with a single esterified fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update glycerol pattern to broadly match a three-carbon chain bearing hydroxy groups
    glycerol_pattern = Chem.MolFromSmarts("[C;R0]([O;R0])[C;R0]([O;R0])[C;R0]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Recognize ester group more exactly as (C(=O)OC) pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;R0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    return True, "Contains glycerol backbone with one esterified fatty acid chain"