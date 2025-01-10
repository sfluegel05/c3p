"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: CHEBI:XXXXX butyrate ester
Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains the butyrate group (CCCC(=O)O-) as part of an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butyrate ester pattern: CCCC(=O)O-[*] (butyrate group as part of an ester)
    butyrate_ester_pattern = Chem.MolFromSmarts("CCCC(=O)O-[*]")
    if not mol.HasSubstructMatch(butyrate_ester_pattern):
        return False, "No butyrate ester group (CCCC(=O)O-[*]) found"

    # Check if the butyrate group is part of an ester (connected to another atom via oxygen)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][*]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Verify that the butyrate group is connected to an ester oxygen
    butyrate_matches = mol.GetSubstructMatches(butyrate_ester_pattern)
    for match in butyrate_matches:
        # Get the oxygen atom in the butyrate group
        oxygen_idx = match[-1]
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        # Check if the oxygen is connected to another atom (ester linkage)
        if oxygen_atom.GetDegree() != 2:
            return False, "Butyrate group not part of an ester linkage"

    return True, "Contains butyrate group (CCCC(=O)O-) as part of an ester"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXX',
        'name': 'butyrate ester',
        'definition': 'Any carboxylic ester where the carboxylic acid component is butyric acid.',
        'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0
}