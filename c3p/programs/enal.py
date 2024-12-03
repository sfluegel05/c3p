"""
Classifies: CHEBI:51688 enal
"""
from rdkit import Chem

def is_enal(smiles: str):
    """
    Determines if a molecule is an enal (alpha,beta-unsaturated aldehyde).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group (C=O at the end of the molecule)
    patt_aldehyde = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
    aldehyde_matches = mol.GetSubstructMatches(patt_aldehyde)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check for conjugated C=C double bond at the alpha,beta position to the aldehyde
    patt_alpha_beta_unsat = Chem.MolFromSmarts('[CX3]=[CX3][CX3H1](=O)')
    alpha_beta_unsat_matches = mol.GetSubstructMatches(patt_alpha_beta_unsat)
    if not alpha_beta_unsat_matches:
        return False, "No alpha,beta-unsaturated aldehyde group found"

    return True, "Molecule is an enal"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51688',
                          'name': 'enal',
                          'definition': 'An alpha,beta-unsaturated aldehyde of '
                                        'general formula R(1)R(2)C=CR(3)-CH=O '
                                        'in which the aldehydic C=O function '
                                        'is conjugated to a C=C double bond at '
                                        'the alpha,beta position.',
                          'parents': ['CHEBI:51718', 'CHEBI:78840']},
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
    'num_false_positives': 1,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 0.9375,
    'recall': 1.0,
    'f1': 0.967741935483871,
    'accuracy': None}