"""
Classifies: CHEBI:38338 aminopyrimidine
"""
from rdkit import Chem

def is_aminopyrimidine(smiles: str):
    """
    Determines if a molecule is an aminopyrimidine (pyrimidine substituted by at least one amino group and its derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aminopyrimidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a pyrimidine ring
    pyrimidine_smarts = 'n1cnccc1'
    pyrimidine = Chem.MolFromSmarts(pyrimidine_smarts)
    if not mol.HasSubstructMatch(pyrimidine):
        return False, "No pyrimidine ring found"

    # Check for amino groups attached to the pyrimidine ring
    amino_smarts = '[NX3;H2,H1;!$(NC=O)]'
    amino_groups = mol.GetSubstructMatches(Chem.MolFromSmarts(amino_smarts))
    
    if not amino_groups:
        return False, "No amino groups found"

    # Check if any amino group is attached to the pyrimidine ring
    pyrimidine_atoms = mol.GetSubstructMatch(pyrimidine)
    pyrimidine_atom_set = set(pyrimidine_atoms)
    
    for amino_group in amino_groups:
        for atom_idx in amino_group:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in pyrimidine_atom_set:
                    return True, "Aminopyrimidine found"

    return False, "No amino groups attached to the pyrimidine ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38338',
                          'name': 'aminopyrimidine',
                          'definition': 'A member of the class of  pyrimidines '
                                        'that is pyrimidine substituted by at '
                                        'least one amino group and its '
                                        'derivatives.',
                          'parents': ['CHEBI:33860', 'CHEBI:39447']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '[00:15:28] Explicit valence for atom # 10 C, 5, is greater than '
             'permitted\n'
             '[00:15:28] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 16,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.8888888888888888,
    'f1': 0.9411764705882353,
    'accuracy': None}