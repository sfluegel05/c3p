"""
Classifies: CHEBI:62732 aromatic ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_ester(smiles: str):
    """
    Determines if a molecule is an aromatic ester.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional group
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester functional group found"

    # Find the ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for match in ester_matches:
        carbonyl_c = match[0]
        ester_o = match[1]

        # Check if the ester oxygen is bonded to an aromatic system
        ester_o_atom = mol.GetAtomWithIdx(ester_o)
        neighbors = ester_o_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() != carbonyl_c and neighbor.GetIsAromatic():
                return True, "Ester linkage is bonded directly to an aromatic system"
    
    return False, "Ester linkage is not bonded directly to an aromatic system"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62732',
                          'name': 'aromatic ester',
                          'definition': 'An ester where the ester linkage is '
                                        'bonded directly to an aromatic '
                                        'system.',
                          'parents': ['CHEBI:33659', 'CHEBI:35701']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 93,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}