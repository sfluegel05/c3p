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

    # Define SMARTS patterns for nucleobases (adenine, guanine, cytosine, thymine, uracil)
    adenine_pattern = Chem.MolFromSmarts('n1c[nH]c2c1ncnc2N')  # Adenine
    guanine_pattern = Chem.MolFromSmarts('O=C1NC(N)=NC2=NC=NC12')  # Guanine
    cytosine_pattern = Chem.MolFromSmarts('O=C1NC=CC(N)=N1')  # Cytosine
    thymine_pattern = Chem.MolFromSmarts('CC1=CN(C)C(=O)NC1=O')  # Thymine
    uracil_pattern  = Chem.MolFromSmarts('O=C1NC=CC(=O)N1')  # Uracil

    nucleobase_patterns = [adenine_pattern, guanine_pattern, cytosine_pattern, thymine_pattern, uracil_pattern]

    # Check for nucleobase
    nucleobase_found = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Define SMARTS patterns for ribose and deoxyribose sugars
    ribose_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O')  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts('[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1')  # Deoxyribose (missing one OH)

    # Check for sugar moiety
    sugar_found = False
    if mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern):
        sugar_found = True
    if not sugar_found:
        return False, "No sugar moiety found"

    # Check for nucleoside linkage (nucleobase attached to sugar via glycosidic bond)
    nucleoside_pattern = Chem.MolFromSmarts('n[cH1][C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O')
    if not mol.HasSubstructMatch(nucleoside_pattern):
        return False, "No nucleoside linkage found"

    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts('OP(=O)(O)[O]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Check if phosphate is attached to sugar hydroxy group
    phosphate_attached = False
    for match in phosphate_matches:
        phosphate_oxygen = mol.GetAtomWithIdx(match[2])  # The oxygen atom connected to the sugar
        for bond in phosphate_oxygen.GetBonds():
            neighbor = bond.GetOtherAtom(phosphate_oxygen)
            if neighbor.GetAtomicNum() == 6:  # Carbon atom in sugar
                phosphate_attached = True
                break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "Phosphate group not attached to sugar"

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
    'attempt': 0,
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