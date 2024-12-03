"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alternating uronic acids and glycosamines
    uronic_acid = False
    glycosamine = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors and 'C' in neighbors and 'O' in [nbr.GetSymbol() for nbr in mol.GetAtomWithIdx(neighbors.index('C')).GetNeighbors()]:
                uronic_acid = True
        if atom.GetSymbol() == 'N' and atom.GetDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'C' in neighbors and 'O' in neighbors:
                glycosamine = True

    if uronic_acid and glycosamine:
        return True, "Molecule contains alternating uronic acids and glycosamines"
    else:
        return False, "Molecule does not contain alternating uronic acids and glycosamines"

    # Check for sulfuric acid esterification
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors:
                return True, "Molecule contains sulfuric acid esterification"

    return False, "Molecule does not meet the criteria for mucopolysaccharide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37395',
                          'name': 'mucopolysaccharide',
                          'definition': 'Any of the group of polysaccharides '
                                        'composed of alternating units from '
                                        'uronic acids and glycosamines, and '
                                        'commonly partially esterified with '
                                        'sulfuric acid.',
                          'parents': ['CHEBI:18085']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 9-10: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}