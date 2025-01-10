"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
# bisbenzylisoquinoline_alkaloid.py
"""
Classifies: bisbenzylisoquinoline alkaloid
Definition: A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units linked by ether bridges. Various structural patterns resulting from additional bridging between the two units by direct carbon-carbon bridging or by methylenedioxy groups are common.
"""

from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges,
    methylenedioxy groups, or direct carbon-carbon bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define isoquinoline core pattern (allows for substitutions)
    isoquinoline_smarts = "c1cc2ccncc2cc1"
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)
    if isoquinoline_pattern is None:
        return False, "Invalid isoquinoline SMARTS pattern"
    
    # Find isoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    num_isoquinolines = len(isoquinoline_matches)
    if num_isoquinolines < 2:
        return False, f"Found {num_isoquinolines} isoquinoline unit(s), need at least 2"
    
    # Assume benzyl group is attached to isoquinoline (flexible matching)
    # Since benzylisoquinoline units are complex, we can consider the isoquinoline unit plus adjacent carbons
    # For simplicity, we check for two isoquinoline units connected via bridges
    
    # Identify connections between isoquinoline units
    isoquinoline_atoms = [set(match) for match in isoquinoline_matches]
    units_connected = False
    for i in range(len(isoquinoline_atoms)):
        for j in range(i+1, len(isoquinoline_atoms)):
            # Check for paths between the two isoquinoline units
            path = Chem.rdmolops.GetShortestPath(mol, list(isoquinoline_atoms[i])[0], list(isoquinoline_atoms[j])[0])
            if path:
                # Check if the path length is reasonable (bridges)
                if len(path) <= 20:  # Adjust as needed for bridging units
                    units_connected = True
                    break
        if units_connected:
            break
    
    if not units_connected:
        return False, "No bridging connections between isoquinoline units found"
    
    return True, "Contains two isoquinoline units connected via bridges"
    
__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'bisbenzylisoquinoline alkaloid',
        'definition': 'A type of benzylisoquinoline alkaloid whose structures are built up of two benzylisoquinoline units linked by ether bridges. Various structural patterns resulting from additional bridging between the two units by direct carbon-carbon bridging or by methylenedioxy groups are common.',
        'parents': ['benzylisoquinoline alkaloid']},
    'config': {   
        'llm_model_name': None,
        'f1_threshold': None,
        'max_attempts': None,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': None,
        'max_negative_in_prompt': None,
        'max_instances_in_prompt': None,
        'test_proportion': None},
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