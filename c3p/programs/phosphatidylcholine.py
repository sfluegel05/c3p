"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the glycerophosphocholine core structure
    glycerophosphocholine_smarts = '[C@@H](COC(=O)[*])(COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)[*]'
    glycerophosphocholine_pattern = Chem.MolFromSmarts(glycerophosphocholine_smarts)
    if not mol.HasSubstructMatch(glycerophosphocholine_pattern):
        return False, "Does not contain glycerophosphocholine core structure"

    # Check for two acyl substituents at positions 1 and 2
    acyl_smarts = 'C(=O)[*]'
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    matches = mol.GetSubstructMatches(acyl_pattern)
    if len(matches) < 2:
        return False, "Does not contain two acyl substituents at positions 1 and 2"

    return True, "Valid phosphatidylcholine structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64482',
                          'name': 'phosphatidylcholine',
                          'definition': 'A glycerophosphocholine that is '
                                        'glycero-3-phosphocholine bearing two '
                                        'acyl substituents at positions 1 and '
                                        '2.',
                          'parents': ['CHEBI:35284', 'CHEBI:36313']},
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
    'num_true_positives': 102,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.9807692307692307,
    'f1': 0.9902912621359222,
    'accuracy': None}