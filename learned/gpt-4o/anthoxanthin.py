"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoids with specific structural characteristics.

    Args:
        smiles (str): SMILES string of molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broad flavonoid core patterns
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("O=c1cc2ccccc2oc1"),       # Flavone pattern
        Chem.MolFromSmarts("Oc1cc2c(O)cc(O)cc2oc1"),  # Flavanol pattern
        Chem.MolFromSmarts("Oc1cc2cc(O)c([O])c(O)c2c(=O)c1"), # Aurone pattern
    ]

    # Check for any flavonoid core structure
    has_core_structure = any(mol.HasSubstructMatch(p) for p in flavonoid_core_patterns)
    if not has_core_structure:
        return False, "Flavonoid core structure not found"

    # Ensure multiple hydroxyl groups, not just one
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups for typical anthoxanthin structure"

    # Check for glycosylation pattern (sugar moieties)
    sugar_patterns = [
        Chem.MolFromSmarts("C1([OX2H][C@H](O)[C@H](O)[C@H](O)C1O)O"),  # Glucose
        Chem.MolFromSmarts("C1O[C@H](C([OX2H])[C@H](O)[C@H](O)[C@H]1O)O"), # Arabinose
    ]
    has_sugar = any(len(mol.GetSubstructMatches(p)) > 0 for p in sugar_patterns)
    
    details = "Flavonoid core with sufficient hydroxyl groups"
    if has_sugar:
        details += " and glycosylation present"
    
    return True, details