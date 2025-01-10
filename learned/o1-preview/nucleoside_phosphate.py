"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleoside (nucleobase + sugar) where one or more of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define general SMARTS patterns for nucleobases
    # Purine base (adenine, guanine)
    purine_pattern = Chem.MolFromSmarts('c1nc[nH]c2c1ncn2')  # General purine ring
    # Pyrimidine base (cytosine, thymine, uracil)
    pyrimidine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1')  # General pyrimidine ring

    # Check for nucleobase (purine or pyrimidine)
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No nucleobase found"

    # Define SMARTS pattern for ribose or deoxyribose sugar
    sugar_pattern = Chem.MolFromSmarts('[C@@H]1O[C@H]([C@@H](O)[C@H](O)[C@H]1O)')  # Ribose ring
    deoxy_sugar_pattern = Chem.MolFromSmarts('[C@@H]1O[C@H]([C@@H](O)[C@H](O)[C@H]1)')  # Deoxyribose ring

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern) and not mol.HasSubstructMatch(deoxy_sugar_pattern):
        return False, "No sugar moiety found"

    # Define SMARTS pattern for nucleoside linkage (nucleobase attached to sugar)
    nucleoside_linkage_pattern = Chem.MolFromSmarts('[nH]1[c,n]%a[c,n]%a[c,n]%a[nH]%a1-[CH]1O[CH][CH][CH][O]1')  # Nucleobase connected to ribose
    if not mol.HasSubstructMatch(nucleoside_linkage_pattern):
        return False, "No nucleoside linkage found"

    # Define SMARTS pattern for phosphate group attached to sugar
    phosphate_pattern = Chem.MolFromSmarts('OP(=O)(O)[O]-[CH2]-')  # Phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group attached to sugar"

    return True, "Contains nucleoside with phosphate(s) attached to sugar"

__metadata__ = {   'chemical_class': {   'id': None,
                          'name': 'nucleoside phosphate',
                          'definition': 'A nucleobase-containing molecular entity that is a nucleoside in which one or more of the sugar hydroxy groups has been converted into a mono- or poly-phosphate. The term includes both nucleotides and non-nucleotide nucleoside phosphates.',
                          'parents': []},
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None}