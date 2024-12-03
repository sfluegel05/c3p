"""
Classifies: CHEBI:60479 lysophosphatidylcholine
"""
from rdkit import Chem

def is_lysophosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a lysophosphatidylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidylcholine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phosphocholine group
    phosphocholine_smarts = '[O-]P(=O)(OCC[N+](C)(C)C)OC'
    phosphocholine = Chem.MolFromSmarts(phosphocholine_smarts)
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Phosphocholine group not found"

    # Check for the glycerol backbone
    glycerol_smarts = '[C@@H](CO)O'
    glycerol = Chem.MolFromSmarts(glycerol_smarts)
    if not mol.HasSubstructMatch(glycerol):
        return False, "Glycerol backbone not found"

    # Check for one acyl group and one hydroxyl group or one hydroxyl group and one acyl group
    acyl_smarts = 'C(=O)O[C@H]'
    acyl = Chem.MolFromSmarts(acyl_smarts)
    hydroxyl_smarts = '[C@H](O)CO'
    hydroxyl = Chem.MolFromSmarts(hydroxyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl)
    
    if len(acyl_matches) + len(hydroxyl_matches) != 2:
        return False, "Incorrect number of acyl or hydroxyl groups"

    return True, "Molecule is a lysophosphatidylcholine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60479',
                          'name': 'lysophosphatidylcholine',
                          'definition': 'An acylglycerophosphocholine '
                                        'resulting from partial hydrolysis of '
                                        'a phosphatidylcholine, which removes '
                                        'one of the fatty acyl groups. The '
                                        'structure is depicted in the image '
                                        'where R1 = acyl, R2 = H or where R1 = '
                                        'H, R2 = acyl.',
                          'parents': ['CHEBI:36313']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 3,
    'num_true_negatives': 15,
    'num_false_negatives': 4,
    'precision': 0.8235294117647058,
    'recall': 0.7777777777777778,
    'f1': 0.7999999999999999,
    'accuracy': None}