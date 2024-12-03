"""
Classifies: CHEBI:22160 acetamides
"""
from rdkit import Chem

def is_acetamides(smiles: str):
    """
    Determines if a molecule is an acetamide (RNHC(=O)CH3).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetamide substructure
    acetamide_smarts = 'NC(=O)C'
    acetamide_substruct = Chem.MolFromSmarts(acetamide_smarts)

    if mol.HasSubstructMatch(acetamide_substruct):
        return True, "Contains acetamide substructure (RNHC(=O)CH3)"
    else:
        return False, "Does not contain acetamide substructure (RNHC(=O)CH3)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22160',
                          'name': 'acetamides',
                          'definition': 'Compounds with the general formula '
                                        'RNHC(=O)CH3.',
                          'parents': ['CHEBI:37622']},
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
    'num_true_positives': 91,
    'num_false_positives': 7,
    'num_true_negatives': 13,
    'num_false_negatives': 2,
    'precision': 0.9285714285714286,
    'recall': 0.978494623655914,
    'f1': 0.9528795811518325,
    'accuracy': None}