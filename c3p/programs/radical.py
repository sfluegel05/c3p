"""
Classifies: CHEBI:26519 radical
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_radical(smiles: str):
    """
    Determines if a molecule is a radical (possessing an unpaired electron).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a radical, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for unpaired electrons
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() > 0:
            return True, f"Atom {atom.GetSymbol()} with index {atom.GetIdx()} has {atom.GetNumRadicalElectrons()} unpaired electron(s)"

    return False, "No unpaired electrons found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26519',
                          'name': 'radical',
                          'definition': 'A molecular entity possessing an '
                                        'unpaired electron.',
                          'parents': ['CHEBI:23367']},
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
    'num_true_positives': 21,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}