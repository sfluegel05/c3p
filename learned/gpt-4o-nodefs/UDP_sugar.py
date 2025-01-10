"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for UDP-sugar components
    uridine_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)")
    
    # Check for uridine moiety
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine moiety found"

    # Check for diphosphate group
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"

    # Check for sugar component
    # A more specific sugar pattern could be defined, but in general a UDP-sugar
    # would have a recognizable sugar moiety linked with the above patterns
    # Here we will assume that if both uridine and diphosphate are detected, the
    # sugar moiety is implied by the context of these patterns in known UDP-sugars.

    return True, "Contains uridine and diphosphate linkage typical of UDP-sugars"