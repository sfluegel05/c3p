"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    N-acylphytosphingosines have a phytosphingosine backbone with an acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a long hydrocarbon chain (>12 carbons, flexible length)
    long_chain_pattern = Chem.MolFromSmarts("C{12,}")  # Using a flexible count-based pattern
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found typical for N-acylphytosphingosine"
    
    # Look for the phytosphingosine backbone (amino alcohol)
    backbone_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](CO)N")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Look for an acyl group bonded to the nitrogen
    acyl_pattern = Chem.MolFromSmarts("N[C;D2](C(=O)[CX4H])")  # Ensuring connection from nitrogen to a carbonyl carbon
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to the nitrogen"
    
    return True, "Contains phytosphingosine backbone with N-acyl group"

# Test the function
example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
result, reason = is_N_acylphytosphingosine(example_smiles)
print(f"Result: {result}, Reason: {reason}")