"""
Classifies: CHEBI:35716 C-nitro compound
"""
from rdkit import Chem

def is_C_nitro_compound(smiles: str):
    """
    Determines if a molecule is a C-nitro compound (a nitro compound with the nitro group attached to a carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a C-nitro compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the nitro group SMARTS pattern
    nitro_pattern = Chem.MolFromSmarts('[CX3](=[NX2+][O-])[N+](=O)[O-]')
    
    if mol.HasSubstructMatch(nitro_pattern):
        return True, "Contains a nitro group attached to a carbon atom"
    else:
        return False, "Does not contain a nitro group attached to a carbon atom"

# Example usage:
# print(is_C_nitro_compound("CC(C)[N+](=O)[O-]")) # Example SMILES for 2-nitropropane
# Output: (True, "Contains a nitro group attached to a carbon atom")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35716',
                          'name': 'C-nitro compound',
                          'definition': 'A nitro compound having the nitro '
                                        'group (-NO2) attached to a carbon '
                                        'atom.',
                          'parents': ['CHEBI:35715', 'CHEBI:72695']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 37-38: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}