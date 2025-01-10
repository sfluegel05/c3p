"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is a furochromene that consists of a furan ring fused with a coumarin.
    The fusion may occur in different ways to give several isomers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for different furanocoumarin cores
    # These patterns capture various fusion isomers between furan and coumarin rings

    patterns = []

    # Linear furanocoumarin (psoralen-type)
    patterns.append(Chem.MolFromSmarts('c1cc2oc(=O)c3ccoc3cc2cc1'))
    
    # Angular furanocoumarin (angelicin-type)
    patterns.append(Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c1ccoc21'))
    
    # Additional patterns to capture other isomers
    patterns.append(Chem.MolFromSmarts('c1cc2oc(=O)c3ccoc3c2c1'))  # Extended linear
    patterns.append(Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c3cocc23'))  # Alternative angular
    patterns.append(Chem.MolFromSmarts('c1cc2c(c1)oc(=O)c3ccoc23'))  # Fusion variation
    patterns.append(Chem.MolFromSmarts('c1ccc2oc(=O)c3ccoc3c2c1'))    # Another variation

    # Check for matches to any of the patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains furanocoumarin core"

    return False, "Does not contain furanocoumarin core"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'furanocoumarin',
        'definition': 'Any furochromene that consists of a furan ring fused with a coumarin. The fusion may occur in different ways to give several isomers.',
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
    }
}