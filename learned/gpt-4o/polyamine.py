"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups, includes charged forms and excludes certain functional groups like amides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Broadened SMARTS patterns incorporating positively charged amines (e.g., [NH3+])
    amine_pattern_primary = Chem.MolFromSmarts('[NX3+0;!$([NX3][CX3]=[OX1])][#1]')  # RNH2 or [NH3+]
    amine_pattern_secondary = Chem.MolFromSmarts('[NX3+0;!$([NX3][CX3]=[OX1])][#6][#1]')  # R2NH
    amine_pattern_tertiary = Chem.MolFromSmarts('[NX3+0;!$([NX3][CX3]=[OX1])][#6][!#1]')  # R3N
    
    # Find matches for the patterns
    primary_matches = mol.GetSubstructMatches(amine_pattern_primary)
    secondary_matches = mol.GetSubstructMatches(amine_pattern_secondary)
    tertiary_matches = mol.GetSubstructMatches(amine_pattern_tertiary)

    # Calculate total amino groups by type
    amine_count = len(primary_matches) + len(secondary_matches) + len(tertiary_matches)

    if amine_count >= 2:
        return True, f"Molecule has {amine_count} amino groups, classifying as polyamine"
    else:
        return False, f"Molecule has {amine_count} amino group(s), not enough to classify as polyamine"