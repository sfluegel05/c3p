"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.

    An amine is defined as 'A compound formally derived from ammonia by replacing
    one, two or three hydrogen atoms by hydrocarbyl groups.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic amine recognition patterns
    basic_amine_patterns = [
        Chem.MolFromSmarts("[NX3;!$(N=*)]"),  # Matches primary/secondary/tertiary amines
        Chem.MolFromSmarts("[NX4+]")          # Matches quaternary ammonium compounds
    ]

    for pattern in basic_amine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Amine group identified"

    return False, "No amine group found"

# Example Usage
smiles_example = "CCC(C)NC"
result, reason = is_amine(smiles_example)
print(f"Result: {result}, Reason: {reason}")

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32988',
                          'name': 'amine',
                          'definition': 'A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.'},
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
    'success': None,
    'best': None,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}