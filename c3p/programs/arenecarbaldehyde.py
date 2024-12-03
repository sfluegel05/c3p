"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde (any aldehyde in which the carbonyl group is attached to an aromatic moiety).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an aldehyde group (formyl group)
    formyl = Chem.MolFromSmarts('[CX3H]=O')
    formyl_matches = mol.GetSubstructMatches(formyl)
    if not formyl_matches:
        return False, "No aldehyde group found"

    # Check for aromaticity
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    if not aromatic_atoms:
        return False, "No aromatic moiety found"

    # Check if the aldehyde group is attached to an aromatic ring
    for match in formyl_matches:
        carbonyl_carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetIdx() in aromatic_atoms:
                return True, "Arenecarbaldehyde detected"

    return False, "Aldehyde group is not attached to an aromatic moiety"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33855',
                          'name': 'arenecarbaldehyde',
                          'definition': 'Any aldehyde in which the carbonyl '
                                        'group is attached to an aromatic '
                                        'moiety.',
                          'parents': ['CHEBI:17478', 'CHEBI:33659']},
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
    'num_true_positives': 16,
    'num_false_positives': 1,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 0.9411764705882353,
    'recall': 1.0,
    'f1': 0.9696969696969697,
    'accuracy': None}