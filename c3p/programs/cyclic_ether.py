"""
Classifies: CHEBI:37407 cyclic ether
"""
from rdkit import Chem

def is_cyclic_ether(smiles: str):
    """
    Determines if a molecule is a cyclic ether.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic ether, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()
    if not rings.AtomRings():
        return False, "No rings found in the molecule"

    # Check if any ring contains an oxygen atom
    for ring in rings.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring):
            return True, "Cyclic ether found"
    
    return False, "No cyclic ether found in the molecule"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37407',
                          'name': 'cyclic ether',
                          'definition': 'Any ether in which the oxygen atom '
                                        'forms part of a ring.',
                          'parents': ['CHEBI:25698', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 1-2: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}