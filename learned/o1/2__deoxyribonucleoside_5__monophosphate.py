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
    # 2'-deoxyribose sugar with correct stereochemistry
    deoxyribose_smarts = """
    [C@H]1([O])O[C@@H](C[C@@H]1O)CO
    """
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)
    if deoxyribose_pattern is None:
        return False, "Error in deoxyribose SMARTS pattern"

    # Phosphate group attached to 5' carbon
    phosphate_smarts = """
    [C@@H]([O])COP(=O)(O)[O]
    """
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if phosphate_pattern is None:
        return False, "Error in phosphate SMARTS pattern"

    # Nucleobases patterns (adenine, guanine, cytosine, thymine)
    adenine_smarts = 'n1(c)ncnc1'
    guanine_smarts = 'n1(c(=O))ncnc1'
    cytosine_smarts = 'n1c(=O)ccn(c1)N'
    thymine_smarts = 'n1c(=O)cc(c1)C'
    nucleobase_smarts_list = [adenine_smarts, guanine_smarts, cytosine_smarts, thymine_smarts]
    nucleobase_patterns = [Chem.MolFromSmarts(s) for s in nucleobase_smarts_list]

    # Check for deoxyribose sugar
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Check for phosphate group at 5' position
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group at 5' position found"

    # Check for nucleobase attached to sugar via N-glycosidic bond
    nucleobase_found = False
    for base_pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(base_pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase attached to the sugar"

    # Combine checks to ensure correct connectivity
    # Full pattern combining sugar, base, and phosphate
    full_pattern_smarts = """
    [$([C@H]1(O)[O][C@@H](C[C@@H]1O)CO)]  # Deoxyribose sugar
    [$([n,c])[c,n]]                        # Nucleobase attached to sugar
    COP(=O)(O)[O]                          # Phosphate group at 5' carbon
    """
    full_pattern = Chem.MolFromSmarts(full_pattern_smarts)
    if full_pattern and mol.HasSubstructMatch(full_pattern):
        return True, "Molecule is a 2'-deoxyribonucleoside 5'-monophosphate"

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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}