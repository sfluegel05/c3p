"""
Classifies: CHEBI:24828 indoles
"""
from rdkit import Chem

def is_indoles(smiles: str):
    """
    Determines if a molecule is an indole (contains an indole skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the indole substructure
    indole_smarts = "c1c[nH]c2ccccc12"
    indole_mol = Chem.MolFromSmarts(indole_smarts)

    if mol.HasSubstructMatch(indole_mol):
        return True, "Contains an indole skeleton"
    
    # Check for other common indole forms
    other_indole_smarts = [
        "c1cnc2ccccc12",  # 1H-indole
        "c1c[nH]c2cccc(c12)",  # 2,3-dihydro-1H-indole
        "c1cnc2cccc(c12)",  # 1H-indole-2,3-dione
        "c1cnc2ccccc2c1=O"  # isatin
    ]
    
    for smarts in other_indole_smarts:
        other_indole_mol = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(other_indole_mol):
            return True, "Contains an indole skeleton (alternative form)"
    
    return False, "Does not contain an indole skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24828',
                          'name': 'indoles',
                          'definition': 'Any compound containing an indole '
                                        'skeleton.',
                          'parents': ['CHEBI:22728']},
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
    'num_true_positives': 181,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 47,
    'precision': 0.9945054945054945,
    'recall': 0.793859649122807,
    'f1': 0.8829268292682927,
    'accuracy': None}