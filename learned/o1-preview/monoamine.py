"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an arylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.

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

    # Define the monoamine pattern: aromatic ring connected via two-carbon chain to amino group
    monoamine_pattern = Chem.MolFromSmarts("a-[CH2]-[CH2]-[NH2]")
    if monoamine_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Check for monoamine substructure
    if mol.HasSubstructMatch(monoamine_pattern):
        return True, "Contains aromatic ring connected via two-carbon chain to amino group"
    else:
        return False, "Does not contain monoamine substructure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'monoamine',
        'definition': 'An arylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain. Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.',
    },
    'message': None,
    'success': True,
    'error': '',
    'stdout': None,
}