"""
Classifies: CHEBI:33238 monoatomic entity
"""
from rdkit import Chem

def is_monoatomic_entity(smiles: str):
    """
    Determines if a molecule is a monoatomic entity (consisting of a single atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoatomic entity, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() == 1:
        return True, "The molecule is a monoatomic entity"
    else:
        return False, "The molecule consists of more than one atom"

# Example usage
# print(is_monoatomic_entity("[Al+3]"))  # Should return (True, "The molecule is a monoatomic entity")
# print(is_monoatomic_entity("CCO"))     # Should return (False, "The molecule consists of more than one atom")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33238',
                          'name': 'monoatomic entity',
                          'definition': 'A monoatomic entity is a molecular '
                                        'entity consisting of a single atom.',
                          'parents': ['CHEBI:33259']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[22:57:06] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[22:57:06] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[22:57:06] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 18,
    'num_false_positives': 0,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}