"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem

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

    # Suppress warnings
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)

    # Define nucleobase patterns (more general, avoiding specific protonation)
    adenine = Chem.MolFromSmarts('n1cnc2ncnc12')  # Adenine
    guanine = Chem.MolFromSmarts('O=C1NC2=NC=NC(N)=C2N1')  # Guanine
    cytosine = Chem.MolFromSmarts('NC1=NC=CC(=O)N1')  # Cytosine
    thymine = Chem.MolFromSmarts('CC1=CN=CN(C1=O)C=O')  # Thymine
    uracil = Chem.MolFromSmarts('O=C1NC=CC(=O)N1')  # Uracil

    nucleobases = [adenine, guanine, cytosine, thymine, uracil]

    # Check for nucleobase
    has_nucleobase = False
    for base in nucleobases:
        if mol.HasSubstructMatch(base):
            has_nucleobase = True
            break
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Define sugar pattern (general five-membered ring with oxygen and hydroxyls)
    sugar_pattern = Chem.MolFromSmarts('C1OC(CO)(CO)C(O)C1')  # Simplified sugar ring

    # Check for sugar moiety connected to nucleobase
    nucleoside_pattern = Chem.MolFromSmarts('*n1cnc2c1ncn2C1OC(O)C(O)C(O)C1')  # Nucleobase connected to sugar
    if not mol.HasSubstructMatch(nucleoside_pattern):
        return False, "No nucleoside linkage (nucleobase attached to sugar) found"

    # Define phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts('OP(=O)(O)O')  # Phosphate group

    # Check for phosphate group(s) attached to sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group attached to sugar"

    # Ensure phosphate is attached to sugar hydroxyl(s)
    phosphate_bonded_to_sugar = False
    sugar_atoms = mol.GetSubstructMatch(sugar_pattern)
    phosphate_atoms = [match[0] for match in phosphate_matches]
    for pa in phosphate_atoms:
        phosphorous = mol.GetAtomWithIdx(pa)
        for neighbor in phosphorous.GetNeighbors():
            if neighbor.GetIdx() in sugar_atoms:
                phosphate_bonded_to_sugar = True
                break
        if phosphate_bonded_to_sugar:
            break
    if not phosphate_bonded_to_sugar:
        return False, "Phosphate not attached to sugar"

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
    'attempt': 2,
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