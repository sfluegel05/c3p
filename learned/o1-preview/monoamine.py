"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an arylalkylamine compound which contains one amino group connected to an aromatic ring by a two-carbon chain.
    Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the monoamine pattern:
    # Aromatic ring connected via two-carbon chain to a nitrogen (primary, secondary, tertiary amine)
    # Exclude nitrogen atoms that are part of a ring
    monoamine_pattern = Chem.MolFromSmarts("""
        [a]                  # Aromatic atom
        -[#6]-[#6]-          # Two-carbon chain
        [NX3;!$(N([#6])[#6]):!$(N=*)]  # Nitrogen with valence 3, not in a ring, not double bonded
    """)
    if monoamine_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Find all matches of the monoamine pattern
    matches = mol.GetSubstructMatches(monoamine_pattern)
    if matches:
        # Check that there is only one amino group connected via two-carbon chain to aromatic ring
        if len(matches) == 1:
            return True, "Contains aromatic ring connected via two-carbon chain to amino group"
        else:
            return False, f"Found {len(matches)} monoamine substructures, expected 1"
    else:
        return False, "Does not contain monoamine substructure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'monoamine',
        'definition': 'An aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain. Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.',
    },
    'message': None,
    'success': True,
    'error': '',
    'stdout': None,
}