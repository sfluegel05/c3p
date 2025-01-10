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
    # Aromatic ring connected via two-carbon chain to a nitrogen (primary, secondary, or tertiary amine)
    # Exclude nitrogen atoms that are part of a ring or double-bonded (imines)
    monoamine_pattern = Chem.MolFromSmarts("[a]-[#6]-[#6]-[NX3;!$(N=C);!R]")

    if monoamine_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Find matches of the monoamine pattern
    matches = mol.GetSubstructMatches(monoamine_pattern)
    if matches:
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