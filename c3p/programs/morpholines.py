"""
Classifies: CHEBI:38785 morpholines
"""
from rdkit import Chem

def is_morpholines(smiles: str):
    """
    Determines if a molecule contains a morpholine ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a morpholine ring, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for morpholine
    morpholine_smiles = "O1CCNCC1"
    morpholine_mol = Chem.MolFromSmiles(morpholine_smiles)

    # Check if the molecule contains the morpholine substructure
    if mol.HasSubstructMatch(morpholine_mol):
        return True, "Molecule contains a morpholine ring"
    else:
        return False, "Molecule does not contain a morpholine ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38785',
                          'name': 'morpholines',
                          'definition': 'Any compound containing morpholine as '
                                        'part of its structure.',
                          'parents': ['CHEBI:46952']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:20:56] SMILES Parse Error: syntax error while parsing: '
             '[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C\\CC\\C(C)=C\x02)C(O)=O\n'
             '[00:20:56] SMILES Parse Error: Failed parsing SMILES '
             "'[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C\\CC\\C(C)=C\x02)C(O)=O' "
             'for input: '
             "'[H][C@@]12OC(=O)C(=C)[C@]1([H])CC\\C(C)=C\\CC\\C(=C\\CC\\C(C)=C\x02)C(O)=O'\n",
    'stdout': '',
    'num_true_positives': 30,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}