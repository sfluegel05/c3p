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
    
    # Identify key patterns: Look for the amino alcohol backbone with a long chain and acyl group
    # Key Pattern 1: Look for a long hydrocarbon chain
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found typical for N-acylphytosphingosine"
    
    # Key Pattern 2: Look for the phytosphingosine backbone (amino alcohol)
    backbone_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](CO)N")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phytosphingosine backbone found"

    # Key Pattern 3: Look for an acyl group bonded to the nitrogen
    acyl_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to the nitrogen"
    
    return True, "Contains phytosphingosine backbone with N-acyl group"

# Test the function
example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)C(O)C(O)C"
result, reason = is_N_acylphytosphingosine(example_smiles)
print(f"Result: {result}, Reason: {reason}")