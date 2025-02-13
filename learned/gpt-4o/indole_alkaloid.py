"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole Alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (a bicyclic structure consisting
    of a benzene ring fused to a pyrrole ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible SMARTS pattern for an indole skeleton
    # Enhance the pattern to accommodate possible substitutions and tautomers
    indole_pattern = Chem.MolFromSmarts('[cH]1ccc2[nH]c([cH]c2c1)')

    # Use alternative SMARTS to capture various indole modifications:
    # A broader pattern matching any nitrogen (n) to account for more substitutions or non-standard tautomeric forms
    indole_pattern_alternative = Chem.MolFromSmarts('c1ccc2[nR]c1ccc2')

    # Match using both primary and alternative SMARTS patterns
    if mol.HasSubstructMatch(indole_pattern) or mol.HasSubstructMatch(indole_pattern_alternative):
        return True, "Contains indole skeleton"
    else:
        return False, "No indole skeleton found"