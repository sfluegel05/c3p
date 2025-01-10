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
        bool: True if the molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate group pattern
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")  # Consider negatively charged sulfate structure
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group found"

    # General sugar-like ring, accounting possible variations in natural sulfolipids
    sugar_ring_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@H](O[C@H]1*)*"),  # generalized sugar like rings
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O[C@@H]1*)*")   # different stereochemistry
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_ring_patterns):
        return False, "No suitable sugar-like ring structure found"

    # Long Chain Fatty Acids: More generic pattern with flexible carbon length
    long_chain_patterns = [
        Chem.MolFromSmarts("C[CH2]{10,}C(=O)O"),                  # Flexible chain length
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")  # Minimum length
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns):
        return False, "No sufficient long fatty acid chain found"

    # If all key features are found
    return True, "Contains a sulfate group, a suitable sugar-like ring structure, and a long fatty acid chain"