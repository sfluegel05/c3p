"""
Classifies: CHEBI:24662 hydroxy-amino acid
"""
from rdkit import Chem

def is_hydroxy_amino_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid core structure: N[C@@H](C)C(O)=O
    amino_acid_core = Chem.MolFromSmarts('N[C@@H](C)C(O)=O')
    if not mol.HasSubstructMatch(amino_acid_core):
        return False, "No amino acid core structure found"

    # Check for hydroxy groups
    hydroxy_group = Chem.MolFromSmarts('[OH]')
    if not mol.HasSubstructMatch(hydroxy_group):
        return False, "No hydroxy groups found"

    return True, "Hydroxy-amino acid structure found"

# Example usage
smiles_examples = [
    "N[C@@H](Cc1cccc(O)c1)C(O)=O",  # L-m-tyrosine
    "OC1CN[C@@H](C1)C(O)=O",  # 4-hydroxy-L-proline
    "N[C@@H](CCCCNO)C(O)=O",  # N(6)-hydroxy-L-lysine
    "NC(CCCN\\C(N)=N/O)C(O)=O",  # N(5)-[(Z)-amino(hydroxyimino)methyl]ornithine
    "NC[C@@H](O)CC[C@H](N)C(O)=O",  # threo-5-hydroxy-L-lysine
    "CC(O)(C[C@H](N)C(O)=O)C(O)=O",  # 4-hydroxy-4-methyl-L-glutamic acid
    "N[C@@H]([C@@H](O)C(O)=O)C(O)=O",  # (3R)-3-hydroxy-L-aspartic acid
    "O=C(O)[C@@H](N(O)O)CCCCCCCSC",  # N,N-dihydroxy-L-pentahomomethionine
    "N[C@H](C(O)=O)c1ccc(O)cc1",  # L-4-hydroxyphenylglycine
    "CSCCCCCCC(N(O)O)C(O)=O",  # N,N-dihydroxytetrahomomethionine
    "N[C@@H]([C@H](O)CC(O)=O)C(O)=O",  # (R)-3-hydroxy-L-glutamic acid
    "N[C@@H](CCCNC(=N)NO)C(O)=O"  # N(5)-[(hydroxyamino)(imino)methyl]-L-ornithine
]

for smiles in smiles_examples:
    result, reason = is_hydroxy_amino_acid(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24662',
                          'name': 'hydroxy-amino acid',
                          'definition': 'A non-proteinogenic alpha-amino acid '
                                        'bearing one or more hydroxy groups at '
                                        'unspecified positions.',
                          'parents': ['CHEBI:83925']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: N[C@@H](Cc1cccc(O)c1)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: OC1CN[C@@H](C1)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: N[C@@H](CCCCNO)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: NC(CCCN\\C(N)=N/O)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: NC[C@@H](O)CC[C@H](N)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: CC(O)(C[C@H](N)C(O)=O)C(O)=O -> True, Hydroxy-amino '
              'acid structure found\n'
              'SMILES: N[C@@H]([C@@H](O)C(O)=O)C(O)=O -> True, Hydroxy-amino '
              'acid structure found\n'
              'SMILES: O=C(O)[C@@H](N(O)O)CCCCCCCSC -> True, Hydroxy-amino '
              'acid structure found\n'
              'SMILES: N[C@H](C(O)=O)c1ccc(O)cc1 -> False, No amino acid core '
              'structure found\n'
              'SMILES: CSCCCCCCC(N(O)O)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n'
              'SMILES: N[C@@H]([C@H](O)CC(O)=O)C(O)=O -> True, Hydroxy-amino '
              'acid structure found\n'
              'SMILES: N[C@@H](CCCNC(=N)NO)C(O)=O -> True, Hydroxy-amino acid '
              'structure found\n',
    'num_true_positives': 11,
    'num_false_positives': 12,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.4782608695652174,
    'recall': 0.9166666666666666,
    'f1': 0.6285714285714286,
    'accuracy': None}