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
    
    # Define the isoalloxazine core with SMARTS
    isoalloxazine_smarts = '[cH]1c(c2nc3c(nc(=O)[nH]c3=O)n(C)c2c(c1)C)=O'
    isoalloxazine_core = Chem.MolFromSmarts(isoalloxazine_smarts)
    if isoalloxazine_core is None:
        return False, "Error in defining isoalloxazine core SMARTS pattern"
    
    # Check for isoalloxazine core
    if not mol.HasSubstructMatch(isoalloxazine_core):
        return False, "No isoalloxazine core found"
    
    # Find matches for the isoalloxazine core
    matches = mol.GetSubstructMatches(isoalloxazine_core)
    if not matches:
        return False, "Isoalloxazine core not matched properly"
    
    # For each match, check for methyl groups at positions 7 and 8
    for match in matches:
        mol_fragment = Chem.PathToSubmol(mol, match)
        
        # Positions of the carbons at positions 7 and 8 in the SMARTS pattern
        pos7_idx = match[1]  # c at index 1 in SMARTS
        pos8_idx = match[5]  # c at index 5 in SMARTS
        
        pos7_atom = mol.GetAtomWithIdx(pos7_idx)
        pos8_atom = mol.GetAtomWithIdx(pos8_idx)
        
        # Check for methyl groups at positions 7 and 8
        methyl7 = False
        methyl8 = False
        
        for neighbor in pos7_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
                methyl7 = True
                break
        for neighbor in pos8_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
                methyl8 = True
                break
        if not methyl7:
            continue  # Try next match
        if not methyl8:
            continue  # Try next match
        
        # Position 10 nitrogen atom in the isoalloxazine core
        pos10_idx = match[7]  # n at index 7 in SMARTS
        pos10_atom = mol.GetAtomWithIdx(pos10_idx)
        
        # Check if position 10 nitrogen has a substituent (non-ring bond)
        substituent = False
        for bond in pos10_atom.GetBonds():
            if not bond.IsInRing():
                substituent = True
                break
        if not substituent:
            continue  # Try next match
        
        # All criteria met
        return True, "Contains dimethylisoalloxazine core with substituent at position 10"
    
    return False, "Molecule does not meet all criteria for a flavin"

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
    'accuracy': None
}