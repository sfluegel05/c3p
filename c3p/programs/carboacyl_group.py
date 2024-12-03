"""
Classifies: CHEBI:37838 carboacyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_carboacyl_group(smiles: str):
    """
    Determines if a molecule is a carboacyl group (formed by loss of at least one OH from the carboxy group of a carboxylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboacyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carbonyl group (C=O) with a wildcard single bond (to represent the loss of OH)
    carboacyl_pattern = Chem.MolFromSmarts('C(=O)[*]')
    if mol.HasSubstructMatch(carboacyl_pattern):
        return True, "Carboacyl group identified"
    
    return False, "Does not match carboacyl group criteria"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37838',
                          'name': 'carboacyl group',
                          'definition': 'A carboacyl group is a group formed '
                                        'by loss of at least one OH from the '
                                        'carboxy group of a carboxylic acid.',
                          'parents': ['CHEBI:22221']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 30,
    'num_false_positives': 15,
    'num_true_negatives': 5,
    'num_false_negatives': 4,
    'precision': 0.6666666666666666,
    'recall': 0.8823529411764706,
    'f1': 0.7594936708860759,
    'accuracy': None}