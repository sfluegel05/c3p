"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphatidyl group
    phosphatidyl_group = Chem.MolFromSmarts('COP(O)(=O)O')
    if not mol.HasSubstructMatch(phosphatidyl_group):
        return False, "No phosphatidyl group found"

    # Check for the presence of the ethanolamine group
    ethanolamine_group = Chem.MolFromSmarts('OCCN')
    if not mol.HasSubstructMatch(ethanolamine_group):
        return False, "No ethanolamine group found"

    # Check if the phosphatidyl group is esterified to the hydroxy group of ethanolamine
    ester_linkage = Chem.MolFromSmarts('COP(O)(=O)OCCN')
    if not mol.HasSubstructMatch(ester_linkage):
        return False, "Phosphatidyl group is not esterified to the hydroxy group of ethanolamine"

    return True, "Valid phosphatidylethanolamine"

# Example usage
smiles = 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COC(=O)CCCCCCCCCCCCCC)(OCCN)(O)=O'
print(is_phosphatidylethanolamine(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16038',
                          'name': 'phosphatidylethanolamine',
                          'definition': 'A class of glycerophospholipids in '
                                        'which a phosphatidyl group is '
                                        'esterified to the hydroxy group of '
                                        'ethanolamine.',
                          'parents': ['CHEBI:36314']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Valid phosphatidylethanolamine')\n",
    'num_true_positives': 114,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}