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
    A monoamine is an aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.
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
    # Aromatic ring connected via two-carbon chain to a nitrogen (amine)
    # Include charged nitrogens ([N;+0,+1])
    # Exclude nitrogens connected to carbonyl groups (not amides)
    # Exclude nitrogen atoms that are part of rings
    monoamine_pattern = Chem.MolFromSmarts("""
    [a]                 # Aromatic atom
    -[#6;!R]            # Non-ring carbon
    -[#6;!R]            # Non-ring carbon
    -[N;X3,X4;+0,+1;!R;!$(N-C=O)]  # Nitrogen not in ring, valence 3 or 4, charge 0 or +1, not connected to C=O
    """)
    if monoamine_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Find matches of the monoamine pattern
    if mol.HasSubstructMatch(monoamine_pattern):
        return True, "Contains aromatic ring connected via two-carbon chain to amino group"
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