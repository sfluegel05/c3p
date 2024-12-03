"""
Classifies: CHEBI:22527 aminopurine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_aminopurine(smiles: str):
    """
    Determines if a molecule is an aminopurine (Any purine having at least one amino substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aminopurine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define purine core SMARTS pattern
    purine_smarts = "c1ncnc2ncnc12"
    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    if not mol.HasSubstructMatch(purine_pattern):
        return False, "No purine core found"
    
    # Define amino group SMARTS pattern
    amino_smarts = "[NX3;H2,H1;!$(NC=O)]"
    amino_pattern = Chem.MolFromSmarts(amino_smarts)
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino substituent found"
    
    return True, "Contains purine core and at least one amino substituent"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22527',
                          'name': 'aminopurine',
                          'definition': 'Any purine having at least one amino '
                                        'substituent.',
                          'parents': ['CHEBI:26401']},
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
    'num_true_positives': 87,
    'num_false_positives': 9,
    'num_true_negatives': 11,
    'num_false_negatives': 3,
    'precision': 0.90625,
    'recall': 0.9666666666666667,
    'f1': 0.9354838709677419,
    'accuracy': None}