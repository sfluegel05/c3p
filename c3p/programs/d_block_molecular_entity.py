"""
Classifies: CHEBI:33676 d-block molecular entity
"""
from rdkit import Chem

def is_d_block_molecular_entity(smiles: str):
    """
    Determines if a molecule is a d-block molecular entity (contains one or more atoms of a d-block element).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a d-block molecular entity, False otherwise
        str: Reason for classification
    """
    d_block_elements = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in d_block_elements:
            return True, f"Contains d-block element: {atom.GetSymbol()}"

    return False, "No d-block elements found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33676',
                          'name': 'd-block molecular entity',
                          'definition': 'A d-block molecular entity is a '
                                        'molecular entity containing one or '
                                        'more atoms of a d-block element.',
                          'parents': ['CHEBI:33497']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:15:21] Explicit valence for atom # 12 O, 3, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 69,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9857142857142858,
    'f1': 0.9928057553956835,
    'accuracy': None}