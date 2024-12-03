"""
Classifies: CHEBI:58958 organosulfate oxoanion
"""
from rdkit import Chem

def is_organosulfate_oxoanion(smiles: str):
    """
    Determines if a molecule is an organosulfate oxoanion (RS(=O)2O(-) where R is an organyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organosulfate oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for the sulfate group pattern
    sulfate_pattern = Chem.MolFromSmarts('OS(=O)(=O)[O-]')
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group found"

    # Check if the sulfate group is attached to an organic group (not inorganic)
    matches = mol.GetSubstructMatches(sulfate_pattern)
    for match in matches:
        sulfate_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in sulfate_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:  # Skip the oxygen atoms of the sulfate group
                if neighbor.GetAtomicNum() in {6, 7, 8, 15, 16}:  # Common atoms in organic groups: C, N, O, P, S
                    return True, "Contains organosulfate oxoanion"
    
    return False, "Sulfate group not attached to an organic group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58958',
                          'name': 'organosulfate oxoanion',
                          'definition': 'An organic anion of general formula '
                                        'RS(=O)2O(-) where R is an organyl '
                                        'group.',
                          'parents': ['CHEBI:25696', 'CHEBI:33482']},
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
    'num_true_positives': 23,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}