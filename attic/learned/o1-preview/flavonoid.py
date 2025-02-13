"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: flavonoid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.

    A flavonoid is defined as any member of the superclass flavonoids whose skeleton is based on 1-benzopyran with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonoid core pattern (1-benzopyran fused ring with aryl at position 2)
    flavonoid_smarts = "c1cc2c(cc1)c(O)cc2-c3cccccc3"  # Simplified pattern for flavonoid core
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)
    if flavonoid_pattern is None:
        return None, "Invalid SMARTS pattern"

    # Check for the flavonoid core pattern
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains flavonoid core structure with aryl substituent at position 2"
    else:
        return False, "Does not contain flavonoid core structure with aryl substituent at position 2"

__metadata__ = {
    'chemical_class': {
        'name': 'flavonoid',
        'definition': 'Any member of the \'superclass\' flavonoids whose skeleton is based on 1-benzopyran with an aryl substituent at position 2. The term was originally restricted to natural products, but is now also used to describe semi-synthetic and fully synthetic compounds.'
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
    'stdout': None
}