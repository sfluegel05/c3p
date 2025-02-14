"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:63517 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside 5'-monophosphate is a 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Deoxyribose sugar pattern (missing OH at 2' position)
    deoxyribose_smarts = '[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O'  # Ribose sugar
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)
    if deoxyribose_pattern is None:
        return False, "Error in deoxyribose SMARTS pattern"

    # Modify the ribose pattern to remove 2'-OH group for deoxyribose
    deoxyribose_no_2oh_smarts = '[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1'  # Remove 2'-OH
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_no_2oh_smarts)
    if deoxyribose_pattern is None:
        return False, "Error in deoxyribose SMARTS pattern"

    # Phosphate group attached to 5' carbon
    phosphate_pattern = Chem.MolFromSmarts('COP(=O)(O)O')
    if phosphate_pattern is None:
        return False, "Error in phosphate SMARTS pattern"

    # Nucleobase attached to sugar
    # General pattern for nucleobase-sugar linkage (N-glycosidic bond)
    nucleobase_pattern = Chem.MolFromSmarts('[$(n1cnc2c1ncnc2),$(c1ccn(c1)n)]')  # Purine or pyrimidine base
    if nucleobase_pattern is None:
        return False, "Error in nucleobase SMARTS pattern"

    # Check for deoxyribose sugar
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Check for phosphate group at 5' position
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group at 5' position found"

    # Check for nucleobase attached to sugar
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase attached to the sugar"

    return True, "Molecule is a 2'-deoxyribonucleoside 5'-monophosphate"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63517',
                          'name': "2'-deoxyribonucleoside 5'-monophosphate",
                          'definition': "A 2'-deoxyribonucleoside monophosphate compound with the phosphate group in the 5'-position.",
                          'parents': ['CHEBI:68499', 'CHEBI:25212']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}