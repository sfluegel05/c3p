"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid with specific structural characteristics.

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

    # Define some SMARTS patterns for core flavonoid structure:
    # Anthoxanthins have a typical flavone core structure
    flavonoid_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)ccc1oc(=O)cc2c1")
    hydroxyl_group_pattern = Chem.MolFromSmarts("[OX2H]")
    methoxy_group_pattern = Chem.MolFromSmarts("O[C;H3]c")
   
    # Check for flavonoid core
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Flavonoid core structure not found"

    # Check for hydroxyl groups - essential for anthoxanthins
    if not mol.HasSubstructMatch(hydroxyl_group_pattern):
        return False, "No hydroxyl groups found, these are essential for anthoxanthins"

    # Methoxy groups are optional; they are common but not required
    methoxy_matches = mol.GetSubstructMatches(methoxy_group_pattern)
    has_methoxy = len(methoxy_matches) > 0

    # Checking elemental composition or large ring structure is beyond SMILES matching
    # Focus on typical anthoxanthin features and report if they include sugar/sulfate moieties

    # Return results based on identified patterns
    if has_methoxy:
        return True, "Flavonoid core with hydroxyl and methoxy groups present"
    else:
        return True, "Flavonoid core with hydroxyl groups present"

    return False, "Does not match the expected features of an anthoxanthin"