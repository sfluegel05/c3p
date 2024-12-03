"""
Classifies: CHEBI:140325 secondary carboxamide
"""
from rdkit import Chem

def is_secondary_carboxamide(smiles: str):
    """
    Determines if a molecule is a secondary carboxamide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary carboxamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxamide group
    carboxamide = False
    primary_amine = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 3 for neighbor in neighbors):
                primary_amine = True
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 1 for neighbor in neighbors):
                carboxamide = True

    if carboxamide and primary_amine:
        return True, "Molecule is a secondary carboxamide"
    else:
        return False, "Molecule is not a secondary carboxamide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140325',
                          'name': 'secondary carboxamide',
                          'definition': 'A carboxamide resulting from the '
                                        'formal condensation of a carboxylic '
                                        'acid with a primary amine; formula '
                                        'RC(=O)NHR(1).',
                          'parents': ['CHEBI:37622']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 2-3: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}