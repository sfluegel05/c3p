"""
Classifies: CHEBI:33554 organosulfonate oxoanion
"""
from rdkit import Chem

def is_organosulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an organosulfonate oxoanion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organosulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonate group (S(=O)(=O)[O-])
    sulfonate_group = Chem.MolFromSmarts('S(=O)(=O)[O-]')
    if not mol.HasSubstructMatch(sulfonate_group):
        return False, "No sulfonate group found"

    # Check if the molecule is an anion
    if not any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms()):
        return False, "Molecule is not an anion"

    # Check if the sulfonate group is part of an organic molecule
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule is not organic"

    return True, "Molecule is an organosulfonate oxoanion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33554',
                          'name': 'organosulfonate oxoanion',
                          'definition': 'An organic anion obtained by '
                                        'deprotonation of the sufonate '
                                        'group(s) of any organosulfonic acid.',
                          'parents': ['CHEBI:25696']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 24-25: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}