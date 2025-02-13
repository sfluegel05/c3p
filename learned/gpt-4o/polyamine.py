"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups (primary, secondary, or tertiary amines).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for primary, secondary, and tertiary amines
    amine_pattern_primary = Chem.MolFromSmarts('[NX3;H2;!$(NC=O)]')  # RNH2
    amine_pattern_secondary = Chem.MolFromSmarts('[NX3;H1;!$(NC=O)]')  # R2NH
    amine_pattern_tertiary = Chem.MolFromSmarts('[NX3;H0;!$(NC=O)]')  # R3N
    
    # Count each type of amines
    amine_matches_primary = mol.GetSubstructMatches(amine_pattern_primary)
    amine_matches_secondary = mol.GetSubstructMatches(amine_pattern_secondary)
    amine_matches_tertiary = mol.GetSubstructMatches(amine_pattern_tertiary)

    # Total amine count
    amine_count = (len(amine_matches_primary) +
                   len(amine_matches_secondary) +
                   len(amine_matches_tertiary))

    # Determine classification
    if amine_count >= 2:
        return True, f"Molecule has {amine_count} amino groups, classifying as polyamine"
    else:
        return False, f"Molecule has {amine_count} amino group(s), not enough to classify as polyamine"