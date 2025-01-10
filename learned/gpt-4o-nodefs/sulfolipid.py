"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate group pattern (corrected)
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O") 
    sulfate_matches = mol.HasSubstructMatch(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate group found"

    # Look for six-membered sugar ring pattern containing oxygens (more specific)
    sugar_ring_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")
    sugar_matches = mol.HasSubstructMatch(sugar_ring_pattern)
    if not sugar_matches:
        return False, "No six-membered ring sugar structure found"

    # Look for long chain fatty acid structures (correct chain length)
    long_chain_pattern = Chem.MolFromSmarts("C[CH2]{12,}[C(=O)O]")  # e.g., longer alkyl chain
    long_chain_matches = mol.HasSubstructMatch(long_chain_pattern)
    if not long_chain_matches:
        return False, "No long fatty acid chain found"

    # If all key features are found
    return True, "Contains a sulfate group, a sugar-like ring structure, and a long fatty acid chain"