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

    # Expand flavonoid core patterns to assess more varieties of benzopyran and other related structures
    flavonoid_core_patterns = [
        # General flavonoid core scaffold (gamma-pyrone structure with two phenyl groups)
        Chem.MolFromSmarts("c1cc2c(c(=O)oc2cc1)-c1ccccc1"),  # Ideal flavone structure for core determination
        Chem.MolFromSmarts("c1cc2oc(=O)c3c(cccc3)cc2c1"),    # Isoflavone-like structure
        Chem.MolFromSmarts("c1cc2c(c1)oc1c3c(c=c(O)cc3)c(=O)oc1c2=O"),  # Chalcone variants
    ]
    
    core_match = any(mol.HasSubstructMatch(core) for core in flavonoid_core_patterns)
    if not core_match:
        return False, "Flavonoid core structure not found"

    # Check for common oxygenation (e.g., multiple hydroxyl or methoxy groups, typical in anthoxanthins)
    oxygenated_patterns = [
        Chem.MolFromSmarts("[O]"),  # General presence of oxygen 
        Chem.MolFromSmarts("[OH]"), # Hydroxyl groups
        Chem.MolFromSmarts("[OX2H]"), # Alcohol or phenol type (e.g., methoxy groups)
    ]
    oxy_matches = sum(mol.HasSubstructMatch(pat) for pat in oxygenated_patterns)
    
    if oxy_matches < 4:
        return False, f"Insufficient oxygen substitutions, found {oxy_matches}. Need at least 4."

    return True, "Contains anthoxanthin characteristics: broad flavonoid scaffold and adequate oxygen substitutions"