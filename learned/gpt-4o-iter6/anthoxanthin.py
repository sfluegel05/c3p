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

    # Broader pattern for flavonoid cores matching variations of benzopyran structures
    flavonoid_core_patterns = [
        # Basic flavonoid core with consideration for flexible structures
        Chem.MolFromSmarts("c1cc2c(c3ccc(O)cc3)c(=O)oc2cc1"),  # benzopyran-4-one core
        Chem.MolFromSmarts("c1cc2c(c1)OC(=O)c3ccccc23"),         # variations of flavone structures
        Chem.MolFromSmarts("c1cc2oc(=O)c3ccccc3cc2c1"),          # isoflavone pattern
    ]
    
    core_match = any(mol.HasSubstructMatch(core) for core in flavonoid_core_patterns)
    if not core_match:
        return False, "Flavonoid core structure not found"

    # Check for sufficient oxygen substitutions (e.g., hydroxyl, methoxy)
    oxygenated_patterns = [Chem.MolFromSmarts("[O]"), Chem.MolFromSmarts("[OH]")]
    oxy_matches = sum(mol.HasSubstructMatch(pat) for pat in oxygenated_patterns)
    
    if oxy_matches < 3:
        return False, f"Insufficient oxygen substitutions, found {oxy_matches}"

    return True, "Contains anthoxanthin characteristics: flavonoid scaffold and adequate oxygen substitutions"