"""
Classifies: CHEBI:24727 hydroxynaphthalene
"""
from rdkit import Chem

def is_hydroxynaphthalene(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthalene (naphthalene carrying one or more hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthalene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for naphthalene core
    naphthalene_smarts = Chem.MolFromSmarts('c1ccc2ccccc2c1')
    if not mol.HasSubstructMatch(naphthalene_smarts):
        return False, "No naphthalene core found"

    # Check for hydroxy groups directly attached to the naphthalene core
    hydroxy_group = Chem.MolFromSmarts('[OH]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_group)
    
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    naphthalene_matches = [match for match in mol.GetSubstructMatches(naphthalene_smarts)]
    for match in hydroxy_matches:
        atom = mol.GetAtomWithIdx(match[0])
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in [idx for submatch in naphthalene_matches for idx in submatch]:
                return True, "Hydroxynaphthalene"

    return False, "No hydroxy groups attached to the naphthalene core found"

# Example usage:
# smiles = "O=C1C2=C(C(O)=CC=C2)C(C=3C=4C(=C(OC)C=CC4)C(O)=CC3)C[C@@H]1C"
# result, reason = is_hydroxynaphthalene(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24727',
                          'name': 'hydroxynaphthalene',
                          'definition': 'Any member of the class of  '
                                        'naphthalenes that is naphthalene '
                                        'carrying one or more hydroxy groups.',
                          'parents': ['CHEBI:25477', 'CHEBI:33853']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 20,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 4,
    'precision': 1.0,
    'recall': 0.8333333333333334,
    'f1': 0.9090909090909091,
    'accuracy': None}