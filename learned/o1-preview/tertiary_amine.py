"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three carbon atoms via single bonds,
    with no hydrogen atoms attached and zero formal charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for tertiary amine nitrogen
    # N with zero formal charge, degree 3, zero hydrogens, connected to three carbons
    tertiary_amine_pattern = Chem.MolFromSmarts("[N;+0;X3;H0;$(N([#6])([#6])[#6])]")

    if tertiary_amine_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Search for matches of the tertiary amine pattern
    matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    if len(matches) == 0:
        return False, "No tertiary amine nitrogen found"

    return True, f"Contains {len(matches)} tertiary amine nitrogen(s)"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32877',
                          'name': 'tertiary amine',
                          'definition': 'A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32874']}}