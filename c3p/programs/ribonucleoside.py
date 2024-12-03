"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (Any nucleoside where the sugar component is D-ribose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-ribose structure as a SMARTS pattern
    ribose_smarts = '[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    
    if ribose_pattern is None:
        return False, "Invalid ribose SMARTS pattern"

    # Check if the molecule contains the D-ribose structure
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No D-ribose component found"

    # Define the nucleoside base structures as SMARTS patterns
    nucleoside_bases = [
        '[nH]1cnc2c1ncnc2',  # Adenine
        'c1cc(=O)[nH]c(=O)n1',  # Uracil
        'c1cc(=O)[nH]c(=O)n1C',  # Thymine
        'c1cc(=O)[nH]cn1',  # Cytosine
        'c1nc(N)nc2c1ncnc2',  # Guanine
    ]
    
    base_patterns = [Chem.MolFromSmarts(base) for base in nucleoside_bases]
    
    if any(base is None for base in base_patterns):
        return False, "Invalid nucleoside base SMARTS pattern(s)"

    # Check if the molecule contains any of the nucleoside base structures
    if not any(mol.HasSubstructMatch(base) for base in base_patterns):
        return False, "No nucleoside base component found"

    return True, "Molecule is a ribonucleoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18254',
                          'name': 'ribonucleoside',
                          'definition': 'Any nucleoside where the sugar '
                                        'component is D-ribose.',
                          'parents': ['CHEBI:33838', 'CHEBI:47019']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}