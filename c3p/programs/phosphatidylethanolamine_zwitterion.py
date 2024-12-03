"""
Classifies: CHEBI:57613 phosphatidylethanolamine zwitterion
"""
from rdkit import Chem

def is_phosphatidylethanolamine_zwitterion(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine zwitterion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine zwitterion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphatidylethanolamine backbone
    phosphate_group = Chem.MolFromSmarts("COP([O-])(=O)OCC[NH3+]")
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphatidylethanolamine backbone found"

    # Ensure the molecule has the zwitterion form
    zwitterion = Chem.MolFromSmarts("COP([O-])(=O)OCC[NH3+]")
    if not mol.HasSubstructMatch(zwitterion):
        return False, "Molecule is not in zwitterion form"

    return True, "Molecule is a phosphatidylethanolamine zwitterion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57613',
                          'name': 'phosphatidylethanolamine zwitterion',
                          'definition': 'The zwitterion of a '
                                        'phosphatidylethanolamine compound '
                                        'formed by proton transfer from the '
                                        'phosphate to the primary amino group.',
                          'parents': ['CHEBI:72823']},
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
    'num_true_positives': 14,
    'num_false_positives': 11,
    'num_true_negatives': 7,
    'num_false_negatives': 4,
    'precision': 0.56,
    'recall': 0.7777777777777778,
    'f1': 0.6511627906976745,
    'accuracy': None}