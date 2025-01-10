"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define more general nucleobase patterns
    # Purine base pattern (adenine, guanine, etc.)
    purine_pattern = Chem.MolFromSmarts('c1ncnc2ncnn12')
    # Pyrimidine base pattern (cytosine, thymine, uracil, etc.)
    pyrimidine_pattern = Chem.MolFromSmarts('c1cncnc1')
    # Include modified bases by allowing substitutions
    nucleobase_pattern = Chem.MolFromSmarts('[#6]1~[#7]~[#6]~[#7]~[#6]1')  # General five-membered ring with Ns and Cs

    # Check for nucleobase presence
    has_nucleobase = False
    nucleobase_matches = mol.GetSubstructMatches(purine_pattern) + mol.GetSubstructMatches(pyrimidine_pattern)
    if nucleobase_matches:
        has_nucleobase = True
    else:
        # Try more general nucleobase pattern
        if mol.HasSubstructMatch(nucleobase_pattern):
            has_nucleobase = True
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Define sugar pattern (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts('C1OC[C@@H](O)[C@@H]1O')
    # Allow for modified sugars by being less specific
    sugar_ring_pattern = Chem.MolFromSmarts('C1OC(C)C(O)C1')  # Five-membered ring with oxygen

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar ring found"

    # Check for nucleoside linkage (nucleobase attached to sugar)
    # Look for an N-C glycosidic bond between nucleobase and sugar
    nucleoside_linkage_pattern = Chem.MolFromSmarts('[nX2]-[CH0]')
    if not mol.HasSubstructMatch(nucleoside_linkage_pattern):
        return False, "No nucleoside linkage (nucleobase attached to sugar) found"

    # Define phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)(O)O')  # Phosphate group

    # Check for phosphate group(s) attached to sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Ensure phosphate is attached to sugar hydroxyl(s)
    phosphate_bonded_to_sugar = False
    sugar_atoms = mol.GetSubstructMatches(sugar_ring_pattern)
    phosphate_atoms = [match[0] for match in phosphate_matches]
    for sugar_match in sugar_atoms:
        for atom_idx in sugar_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 15:  # Phosphorus
                    phosphate_bonded_to_sugar = True
                    break
            if phosphate_bonded_to_sugar:
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
    'attempt': 3,
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