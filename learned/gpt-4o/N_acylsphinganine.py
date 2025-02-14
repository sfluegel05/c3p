"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is defined as a ceramide consisting of sphinganine
    in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sphinganine backbone: C(N)(CO)[C@H](O)CCCC
    sphinganine_pattern = Chem.MolFromSmarts("N[C@H](CO)[C@H](O)CCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # N-acyl linkage: N-C(=O)-C (long chain) - here we make it flexible to include long chains
    n_acyl_pattern = Chem.MolFromSmarts("N[C@H](CO)[C@H](O)[C;R0]-[C;R0](=O)")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl substitution found"

    # Check for long chains with or without hydroxyl groups attached to the carbonyl
    long_chain_patterns = [
        Chem.MolFromSmarts("[C;R0](=O)C([OH])CCCCCCCCCCCCCCCC"),
        Chem.MolFromSmarts("[C;R0](=O)C([OH])CCCCCCCCCCCCCCCCC"),
        Chem.MolFromSmarts("[C;R0](=O)CCCCCCCCCCCCCCCCCCC"),  # Longer carbon chain without hydroxyl
        Chem.MolFromSmarts("[C;R0](=O)CCCCCCCCCCCCCCCC")      # Shorter chain with less termination
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns):
        return False, "No suitable N-acyl long fatty acyl chain detected"

    return True, "Contains sphinganine backbone with N-acyl substitution"