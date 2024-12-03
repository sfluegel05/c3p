"""
Classifies: CHEBI:33560 p-block element atom
"""
from rdkit import Chem

def is_p_block_element_atom(smiles: str):
    """
    Determines if a molecule is a p-block element atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a p-block element atom, False otherwise
        str: Reason for classification
    """
    p_block_elements = {'B', 'C', 'N', 'O', 'F', 'Al', 'Si', 'P', 'S', 'Cl', 'Ga', 'Ge', 'As', 'Se', 'Br', 'In', 'Sn', 'Sb', 'Te', 'I', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    atoms = mol.GetAtoms()
    if len(atoms) != 1:
        return False, "Not a single atom"

    atom = atoms[0]
    symbol = atom.GetSymbol()

    if symbol in p_block_elements:
        return True, f"{symbol} is a p-block element atom"
    else:
        return False, f"{symbol} is not a p-block element atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33560',
                          'name': 'p-block element atom',
                          'definition': 'Any main group element atom belonging '
                                        'to the p-block of the periodic table.',
                          'parents': ['CHEBI:33318']},
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
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9166666666666666,
    'f1': 0.9565217391304348,
    'accuracy': None}