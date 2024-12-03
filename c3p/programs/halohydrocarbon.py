"""
Classifies: CHEBI:24472 halohydrocarbon
"""
from rdkit import Chem

def is_halohydrocarbon(smiles: str):
    """
    Determines if a molecule is a halohydrocarbon (a compound derived from a hydrocarbon by replacing a hydrogen atom with a halogen atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a halohydrocarbon, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    halogens = {'F', 'Cl', 'Br', 'I'}
    has_carbon = False
    has_halogen = False
    has_hydrogen = False
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            has_carbon = True
        if atom.GetSymbol() in halogens:
            has_halogen = True
        if atom.GetSymbol() == 'H':
            has_hydrogen = True
    
    if has_carbon and has_halogen:
        return True, "Molecule is a halohydrocarbon"
    elif not has_carbon:
        return False, "Molecule does not contain a carbon atom"
    elif not has_halogen:
        return False, "Molecule does not contain a halogen atom"
    
    return False, "Molecule does not meet halohydrocarbon criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24472',
                          'name': 'halohydrocarbon',
                          'definition': 'A compound derived from a hydrocarbon '
                                        'by replacing a hydrogen atom with a '
                                        'halogen atom.',
                          'parents': ['CHEBI:17792']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 15,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}