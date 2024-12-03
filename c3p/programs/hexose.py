"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose (six-carbon monosaccharide with an aldehyde group at position 1 or a ketone group at position 2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for six carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 6:
        return False, "Molecule does not have six carbon atoms"

    # Check for aldehyde group at position 1 or ketone group at position 2
    has_aldehyde = False
    has_ketone = False
    aldehyde_pattern = Chem.MolFromSmarts('[CH1](=O)')
    ketone_pattern = Chem.MolFromSmarts('[CH2](=O)[CH2]')
    
    if mol.HasSubstructMatch(aldehyde_pattern):
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            if mol.GetAtomWithIdx(match[0]).GetIdx() == 0:
                has_aldehyde = True
                break

    if mol.HasSubstructMatch(ketone_pattern):
        for match in mol.GetSubstructMatches(ketone_pattern):
            if mol.GetAtomWithIdx(match[0]).GetIdx() == 1:
                has_ketone = True
                break

    if has_aldehyde:
        return True, "Molecule is an aldohexose"
    elif has_ketone:
        return True, "Molecule is a ketohexose"
    else:
        return False, "Molecule does not have an aldehyde group at position 1 or a ketone group at position 2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18133',
                          'name': 'hexose',
                          'definition': 'Any six-carbon monosaccharide which '
                                        'in its linear form contains either an '
                                        'aldehyde group at position 1 '
                                        '(aldohexose) or a ketone group at '
                                        'position 2 (ketohexose).',
                          'parents': ['CHEBI:35381']},
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
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 30,
    'precision': 1.0,
    'recall': 0.03225806451612903,
    'f1': 0.0625,
    'accuracy': None}