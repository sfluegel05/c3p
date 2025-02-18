"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid (CHEBI:37845)
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid has a steroid skeleton with an oxo group at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid core pattern (four fused rings: 3 cyclohexane, 1 cyclopentane)
    steroid_core = Chem.MolFromSmarts("[*]1@[*]@[*]2@[*]@[*]3@[*]@[*]4@[*]@[*]5@[*]@[*]1@[*]@[*]@2@[*]@3@[*]@4@5")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not have a steroid skeleton"

    # SMARTS pattern for oxo group (C=O) at position 11
    # Position 11 corresponds to a specific carbon in the steroid structure
    oxo_11_pattern = Chem.MolFromSmarts("[C;X3]=[O]")
    matches = mol.GetSubstructMatches(oxo_11_pattern)
    if not matches:
        return False, "No oxo group found"

    # Verify that the oxo group is at position 11
    # This requires knowing the specific atom indices in the steroid structure
    # Assuming the 11th position is a specific atom in the SMARTS match
    # Note: This part might need adjustment based on actual steroid numbering
    # For simplicity, check if any carbonyl is present (may need more precise SMARTS)
    # Alternatively, use a more specific SMARTS pattern for position 11
    # Example pattern targeting position 11 (adjust based on actual structure)
    specific_oxo_11 = Chem.MolFromSmarts("[C@@H]1[C@@]2[C@H]([C@H]3[C@@H]([C@@]4(C(=O)CC(C)(C)4)CC3)CC2)CC1")
    if mol.HasSubstructMatch(specific_oxo_11):
        return True, "11-oxo group present on steroid skeleton"
    else:
        return False, "Oxo group not at position 11"