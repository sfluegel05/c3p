"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion (conjugate base of a fatty acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group [C(=O)[O-]]
    carboxylate_group = Chem.MolFromSmarts("[C](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_group):
        return False, "No carboxylate group found"

    # Check if the carboxylate group is attached to a long aliphatic chain (at least 4 carbons)
    aliphatic_chain = Chem.MolFromSmarts("[C;R0;!$(C=*)]!@[C;R0;!$(C=*)]!@[C;R0;!$(C=*)]!@[C;R0;!$(C=*)]")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No long aliphatic chain found"

    return True, "Molecule is a fatty acid anion"

# Example usage:
smiles = "CC(C)[C@@H](O)CC[C@@H](C)CC([O-])=O"
print(is_fatty_acid_anion(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28868',
                          'name': 'fatty acid anion',
                          'definition': 'The conjugate base of a fatty acid, '
                                        'arising from deprotonation of the '
                                        'carboxylic acid group of the '
                                        'corresponding fatty acid.',
                          'parents': ['CHEBI:18059', 'CHEBI:35757']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is a fatty acid anion')\n",
    'num_true_positives': 73,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 27,
    'precision': 0.9864864864864865,
    'recall': 0.73,
    'f1': 0.839080459770115,
    'accuracy': None}