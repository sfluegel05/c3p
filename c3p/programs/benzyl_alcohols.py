"""
Classifies: CHEBI:22743 benzyl alcohols
"""
from rdkit import Chem

def is_benzyl_alcohols(smiles: str):
    """
    Determines if a molecule is a benzyl alcohol (contains a phenylmethanol skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzyl alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the phenylmethanol substructure
    phenylmethanol_smiles = "c1ccccc1CO"
    phenylmethanol = Chem.MolFromSmiles(phenylmethanol_smiles)

    if mol.HasSubstructMatch(phenylmethanol):
        return True, "Molecule contains a phenylmethanol skeleton"
    else:
        return False, "Molecule does not contain a phenylmethanol skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22743',
                          'name': 'benzyl alcohols',
                          'definition': 'Compounds containing a phenylmethanol '
                                        'skeleton.',
                          'parents': ['CHEBI:33854']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 4,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7333333333333333,
    'recall': 1.0,
    'f1': 0.846153846153846,
    'accuracy': None}