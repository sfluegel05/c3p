"""
Classifies: CHEBI:24780 imidazoles
"""
from rdkit import Chem

def is_imidazoles(smiles: str):
    """
    Determines if a molecule is an imidazole or its derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an imidazole or its derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the imidazole core
    imidazole_core = Chem.MolFromSmarts('n1cncc1')
    if imidazole_core is None:
        return False, "Failed to create imidazole core pattern"

    # Check if the molecule contains the imidazole core
    if not mol.HasSubstructMatch(imidazole_core):
        return False, "No imidazole core found"

    return True, "Imidazole or its derivative"

# Example usage
smiles = "N1(C=C(N=C1)C#N)C"  # 1-Methyl-1H-imidazole-4-carbonitrile
result, reason = is_imidazoles(smiles)
print(result, reason)  # Expected output: True, "Imidazole or its derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24780',
                          'name': 'imidazoles',
                          'definition': 'A five-membered organic heterocycle '
                                        'containing two nitrogen atoms at '
                                        'positions 1 and 3, or any of its '
                                        'derivatives; compounds containing an '
                                        'imidazole skeleton.',
                          'parents': ['CHEBI:23677']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Imidazole or its derivative\n',
    'num_true_positives': 44,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 4,
    'precision': 1.0,
    'recall': 0.9166666666666666,
    'f1': 0.9565217391304348,
    'accuracy': None}