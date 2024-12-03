"""
Classifies: CHEBI:25716 organothiophosphorus compound
"""
from rdkit import Chem

def is_organothiophosphorus_compound(smiles: str):
    """
    Determines if a molecule is an organothiophosphorus compound (contains a phosphorus-sulfur bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organothiophosphorus compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of phosphorus-sulfur (P-S) bond
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetSymbol() == 'P' and atom2.GetSymbol() == 'S') or (atom1.GetSymbol() == 'S' and atom2.GetSymbol() == 'P'):
            return True, "Contains phosphorus-sulfur (P-S) bond"

    return False, "Does not contain phosphorus-sulfur (P-S) bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25716',
                          'name': 'organothiophosphorus compound',
                          'definition': 'An organothiophosphorus compound is '
                                        'an organophosphorus compound which '
                                        'contains a phosphorus-sulfur bond.',
                          'parents': ['CHEBI:25710', 'CHEBI:26835']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.7692307692307693,
    'f1': 0.8695652173913044,
    'accuracy': None}