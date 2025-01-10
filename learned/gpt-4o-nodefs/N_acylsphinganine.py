"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.

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

    # Redefine sphinganine backbone pattern with broader matching
    sphinganine_backbone = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)NC(=O)")
    if not mol.HasSubstructMatch(sphinganine_backbone):
        return False, "No sphinganine backbone found"

    # Redefine acyl linkage pattern to allow longer aliphatic chains
    acyl_chain_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[C,CX4]*")  # Allow flexible carbon counts
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No N-acyl linkage found"

    # Check for optional headgroup, if present (sugar moieties or similar)
    possible_headgroup = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)C[C@H]1O")  # Simplified sugar example
    if mol.HasSubstructMatch(possible_headgroup):
        return True, "Contains sphinganine backbone with N-acyl linkage and possible sugar headgroup"

    return True, "Contains sphinganine backbone with N-acyl linkage"