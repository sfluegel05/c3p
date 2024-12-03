"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-lactam ring SMARTS pattern
    beta_lactam_smarts = "[C;R1]1=[O;R1][N;R1][C;R1]1"
    beta_lactam = Chem.MolFromSmarts(beta_lactam_smarts)

    if mol.HasSubstructMatch(beta_lactam):
        return True, "Contains a beta-lactam ring"
    else:
        return False, "Does not contain a beta-lactam ring"

# Example usage
smiles = "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@H](C3=CC=C(C=C3)O)N)C(=O)O)C"
result, reason = is_beta_lactam_antibiotic(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27933',
                          'name': 'beta-lactam antibiotic',
                          'definition': 'An organonitrogen heterocyclic '
                                        'antibiotic that contains a '
                                        'beta-lactam ring.',
                          'parents': ['CHEBI:25558', 'CHEBI:35627']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'False Does not contain a beta-lactam ring\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 32,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}