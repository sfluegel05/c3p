"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:35277 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is an organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary ammonium ion
    # Nitrogen with positive charge, two hydrogens, bonded to two carbons
    pattern = Chem.MolFromSmarts("[N+H2]([C])[C]")
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No secondary ammonium ion group found"

    # Check each match to ensure the nitrogen atom has only two carbon attachments
    for match in matches:
        n_idx = match[0]  # Index of the nitrogen atom
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Check that the nitrogen has a +1 formal charge
        if n_atom.GetFormalCharge() != 1:
            continue
        # Check that the nitrogen is bonded to exactly two carbons
        neighbor_carbons = [atom.GetAtomicNum() for atom in n_atom.GetNeighbors() if atom.GetAtomicNum() == 6]
        if len(neighbor_carbons) != 2:
            continue
        return True, "Contains secondary ammonium ion group ([N+H2]([C])[C])"

    return False, "No secondary ammonium ion group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35277',
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['CHEBI:35274']
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
    'accuracy': None
}