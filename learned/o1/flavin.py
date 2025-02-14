"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione)
    skeleton, with a substituent on the 10 position.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the isoalloxazine core with atom maps for positions 7, 8, and 10
    isoalloxazine_smarts = """
    [#6][$([#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)]
    -c2ncnc3nc(=O)[nH]c(=O)c23
    """
    isoalloxazine_core = Chem.MolFromSmarts(isoalloxazine_smarts)
    if isoalloxazine_core is None:
        return None, "Error in defining isoalloxazine core SMARTS pattern"
    
    # Check for isoalloxazine core
    if not mol.HasSubstructMatch(isoalloxazine_core):
        return False, "No isoalloxazine core found"
    
    # Get the matching substructure
    matches = mol.GetSubstructMatch(isoalloxazine_core)
    if not matches:
        return False, "No match for the isoalloxazine core"
    
    # Map indices to identify positions
    match_atom_indices = {atom.GetIdx() for atom in mol.GetSubstructMatch(isoalloxazine_core)}
    atom_positions = {}
    for atom in isoalloxazine_core.GetAtoms():
        atom_map_num = atom.GetAtomMapNum()
        if atom_map_num:
            atom_positions[atom_map_num] = matches[atom.GetIdx()]
    
    # Check for methyl groups at positions 7 and 8
    position7_idx = matches[0]  # Assuming position 7 is the first atom in pattern
    position8_idx = matches[1]  # Assuming position 8 is the second atom in pattern
    position7_atom = mol.GetAtomWithIdx(position7_idx)
    position8_atom = mol.GetAtomWithIdx(position8_idx)
    
    # Check if atoms at positions 7 and 8 are carbons
    if position7_atom.GetAtomicNum() != 6 or position8_atom.GetAtomicNum() != 6:
        return False, "Positions 7 and 8 are not carbon atoms"
    
    # Check for methyl groups at positions 7 and 8
    methyl7 = False
    methyl8 = False
    for neighbor in position7_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
            methyl7 = True
            break
    for neighbor in position8_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
            methyl8 = True
            break
    if not (methyl7 and methyl8):
        return False, "Methyl groups at positions 7 and 8 not found"
    
    # Check for substituent at position 10 (nitrogen atom in central ring)
    position10_atom = None
    for atom in isoalloxazine_core.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetDegree() == 3:
            # Nitrogen with three bonds is likely position 10
            position10_idx = matches[atom.GetIdx()]
            position10_atom = mol.GetAtomWithIdx(position10_idx)
            break
    if position10_atom is None:
        return False, "Nitrogen at position 10 not found"
    
    # Check if nitrogen at position 10 has a substituent (excluding ring bonds)
    substituents = 0
    for bond in position10_atom.GetBonds():
        if not bond.IsInRing():
            substituents += 1
    if substituents == 0:
        return False, "No substituent found at position 10"
    
    return True, "Contains dimethylisoalloxazine core with substituent at position 10"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'flavin',
        'definition': 'A derivative of the dimethylisoalloxazine (7,8-dimethylbenzo[g]pteridine-2,4(3H,10H)-dione) skeleton, with a substituent on the 10 position.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
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
    'accuracy': None
}