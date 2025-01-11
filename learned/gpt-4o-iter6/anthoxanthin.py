"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments characterized by a benzopyran-4-one structure often with oxygen substitutions.

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

    # Flavonoid basic skeleton: benzopyran and phenyl ring
    # Searching benzopyran-4-one (flavone) core structure with variable substitution
    flavonoid_core_pattern = Chem.MolFromSmarts("c1cc(c2c(c1)C(=O)c3c(occ3)c2)")
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "Flavonoid core structure (benzopyran-4-one) not found"
    
    # Check for common oxygen substitutions (like OH and OCH3)
    oxy_substitutions = Chem.MolFromSmarts("[OH0,O]")
    oxy_matches = len(mol.GetSubstructMatches(oxy_substitutions))
    if oxy_matches < 2:
        return False, f"Insufficient oxygen substitutions, found {oxy_matches}"

    # If necessary, further checks can consider specific methoxy groups, hydrophilic groups
    # that enhance water solubility, which are common in anthoxanthins.
    
    return True, "Contains anthoxanthin characteristics: flavonoid scaffold and oxygen substitutions"