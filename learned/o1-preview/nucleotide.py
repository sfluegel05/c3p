"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33504 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substructures for nucleobases
    purine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')
    pyrimidine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1') # uracil
    cytosine_pattern = Chem.MolFromSmarts('c1cc(nc(N)n1)=O')
    thymine_pattern = Chem.MolFromSmarts('c1c(N)cc(=O)[nH]c1=O')
    hypoxanthine_pattern = Chem.MolFromSmarts('n1c2ncnc2c(=O)[nH]c1')
    xanthine_pattern = Chem.MolFromSmarts('n1c2ncnc2c(=O)[nH]c(=O)c1')

    # Combine nucleobase patterns
    nucleobase_patterns = [purine_pattern, pyrimidine_pattern, cytosine_pattern,
                           thymine_pattern, hypoxanthine_pattern, xanthine_pattern]

    has_nucleobase = any(mol.HasSubstructMatch(pat) for pat in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Define pattern for sugar (allowing for ribose and deoxyribose)
    sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O1')  # Ribose
    deoxy_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](CO)[C@@H](CO)O1')  # Deoxyribose
    modified_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO*)O1')  # Allows modifications

    # Check for sugar
    has_sugar = mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxy_sugar_pattern) or mol.HasSubstructMatch(modified_sugar_pattern)
    if not has_sugar:
        return False, "No pentose sugar found"

    # Check for nucleoside linkage (nucleobase-sugar bond)
    nucleoside_bond_pattern = Chem.MolFromSmarts('[$([nH]),$(n),$(N)]1[c,n]n[c,n][c,n][c,n]1[C@H]2O[C@H]([C@@H](O)[C@H](O)[C@@H]2O)')  # N-glycosidic bond
    has_nucleoside_linkage = mol.HasSubstructMatch(nucleoside_bond_pattern)
    if not has_nucleoside_linkage:
        return False, "No nucleoside linkage found"

    # Check for phosphate group attached to the sugar's 3' or 5' hydroxyl group
    phosphate_pattern = Chem.MolFromSmarts('O[P](=O)(O)[O][C@H]')  # Phosphate ester linkage
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group attached to the sugar's 3' or 5' hydroxyl group"

    return True, "Molecule is a nucleotide with nucleobase, sugar, and phosphate group"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33504',
                              'name': 'nucleotide',
                              'definition': 'A nucleoside phosphate resulting from the condensation of the 3\' or 5\' hydroxy group of a nucleoside with phosphoric acid.',
                              'parents': ['CHEBI:26161']},
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