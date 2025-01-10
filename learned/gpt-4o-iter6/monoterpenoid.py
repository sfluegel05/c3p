"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid generally inherits features from a C10 monoterpene backbone, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Allow some flexibility in carbon count due to rearrangement/modification
    if not (8 <= c_count <= 12):
        return False, f"Carbon count {c_count} not typical for monoterpenoids (expected ~10)"

    # Look for typical monoterpenoid motifs (simplified rule example)
    # This is a placeholder for a more informed pattern or set of checks
    # Look for a simple cyclohexene pattern indicative of a monoterpene origin:
    cyclohexene_pattern = Chem.MolFromSmarts("C1=CCCCC1")
    if not mol.HasSubstructMatch(cyclohexene_pattern):
        return False, "No typical monoterpenoid backbone detected"

    # Presence of terpenoid functional groups, example alcohol, ketone or ester
    # Placeholder patterns for specificity can be defined
    functional_groups = ["[CX3](=O)", "[OX2H]", "[CX3](=O)O[CX3]"]  # ketone, alcohol, ester
    for fg_smarts in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if mol.HasSubstructMatch(fg_pattern):
            return True, "Contains structural features typical of monoterpenoids"

    return False, "Does not have identifiable monoterpenoid features"