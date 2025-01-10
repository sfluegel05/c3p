"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a variety of core structures and multiple oxygen substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expand flavonoid core patterns to cover varied anthoxanthin core structures
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("c1cc2oc(=O)cc(c2c1)-c1ccccc1"),  # General flavone structure
        Chem.MolFromSmarts("c1cc2oc(=O)cc(c2c1)-c1c([!#1])c([!#1])c([!#1])c([!#1])c1"),  # Extended core
        Chem.MolFromSmarts("c1cc2nc(=O)cc(c2c1)-c1ccccc1O"),  # Other chromone-like derivatives
    ]
    
    core_match = any(mol.HasSubstructMatch(core) for core in flavonoid_core_patterns)
    if not core_match:
        return False, "Flavonoid core structure not found"

    # Check for common high degree of oxygenation
    oxygenated_patterns = [
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl groups
        Chem.MolFromSmarts("[OX2H]"),  # Methoxy groups
        Chem.MolFromSmarts("c-O"),  # General oxygen bonded to aromatic carbon
    ]
    oxy_matches = sum(mol.HasSubstructMatch(pat) for pat in oxygenated_patterns)
    
    if oxy_matches < 4:
        return False, f"Insufficient oxygen substitutions, found {oxy_matches}. Need at least 4."

    return True, "Contains anthoxanthin characteristics with flavonoid scaffold and adequate oxygen substitutions"