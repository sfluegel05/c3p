"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:33853 aromatic primary alcohol
Definition: Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a primary alcohol group (-CH2OH) attached to a carbon in an aromatic ring or directly bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the primary alcohol pattern (-CH2OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group (-CH2OH) found"

    # Find all matches for the primary alcohol group
    primary_alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    # Check if any of the primary alcohol groups are attached to an aromatic carbon or a carbon directly bonded to an aromatic ring
    for match in primary_alcohol_matches:
        carbon_idx = match[0]  # Index of the carbon in the -CH2OH group
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Check if the carbon is part of an aromatic ring
        if carbon_atom.GetIsAromatic():
            return True, "Primary alcohol group (-CH2OH) attached to an aromatic carbon"
        
        # Check if the carbon is directly bonded to an aromatic ring
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                # Ensure the aromatic ring is directly bonded to the carbon bearing the -CH2OH group
                return True, "Primary alcohol group (-CH2OH) attached to a carbon directly bonded to an aromatic ring"

    return False, "Primary alcohol group (-CH2OH) not attached to an aromatic carbon or a carbon directly bonded to an aromatic ring"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33853',
        'name': 'aromatic primary alcohol',
        'definition': 'Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.',
        'parents': ['CHEBI:33852', 'CHEBI:25805']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}