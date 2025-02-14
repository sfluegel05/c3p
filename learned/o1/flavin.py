"""
Classifies: CHEBI:30527 flavin
"""
"""
Classifies: flavin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Define the isoalloxazine core without substitutions
    # This pattern matches the tricyclic ring system of isoalloxazine
    isoalloxazine_smarts = 'O=C1NC2=NC(=O)N=C2N(C3=CC=CC=C13)'
    isoalloxazine_core = Chem.MolFromSmarts(isoalloxazine_smarts)
    if isoalloxazine_core is None:
        return None, "Error in defining isoalloxazine core SMARTS pattern"
    
    # Check for isoalloxazine core
    if not mol.HasSubstructMatch(isoalloxazine_core):
        return False, "No isoalloxazine core found"
    
    # Get the matching substructure
    matches = mol.GetSubstructMatches(isoalloxazine_core)
    if not matches:
        return False, "No match for the isoalloxazine core"
    
    # Verify methyl groups at positions 7 and 8
    # Positions 7 and 8 correspond to the carbons adjacent to the benzene ring
    methyl_positions_smarts = 'Cc1cc2nc3c(nc(=O)[nH]c3=O)n(c2cc1C)[N]'
    methyl_pattern = Chem.MolFromSmarts(methyl_positions_smarts)
    if methyl_pattern is None:
        return None, "Error in defining methyl positions SMARTS pattern"
    
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Methyl groups at positions 7 and 8 not found"
    
    # Check for substituent at position 10
    # Position 10 is the nitrogen atom in the pyrimidine ring
    # Find the nitrogen atom at position 10
    nitrogen10_smarts = '[nH]c1ncnc2c1n(C)c(=O)[nH]c2=O'
    nitrogen10_pattern = Chem.MolFromSmarts(nitrogen10_smarts)
    if nitrogen10_pattern is None:
        return None, "Error in defining position 10 nitrogen SMARTS pattern"
    
    match = mol.GetSubstructMatch(nitrogen10_pattern)
    if not match:
        return False, "Nitrogen at position 10 not found"
    
    # Get the atom index of nitrogen at position 10
    nitrogen_idx = match[0]
    nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
    
    # Check if nitrogen has a substituent (connected to an atom outside the core)
    substituent_found = False
    for neighbor in nitrogen_atom.GetNeighbors():
        if neighbor.GetIdx() not in match:
            substituent_found = True
            break
    
    if not substituent_found:
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
    'accuracy': None
}