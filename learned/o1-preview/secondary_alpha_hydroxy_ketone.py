"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone contains a ketone group (C=O) and a hydroxyl group (OH)
    attached to the alpha carbon, which is secondary (attached to one hydrogen, one carbonyl carbon,
    one hydroxyl group, and one other carbon group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary alpha-hydroxy ketone
    # Pattern explanation:
    # [#6](=O)          - Carbonyl carbon (C=O) attached to...
    # [#6;H1](-[OH])    - Alpha carbon [#6;H1] with one hydrogen and a hydroxyl group (-OH)
    # -[#6]             - Alpha carbon connected to another carbon (organyl group)
    pattern = Chem.MolFromSmarts("[#6](=O)-[#6;H1](-[OH])-[#6]")

    # Check for substructure match
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule matches secondary alpha-hydroxy ketone pattern"
    else:
        return False, "Molecule does not match secondary alpha-hydroxy ketone pattern"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'secondary alpha-hydroxy ketone',
        'definition': 'An alpha-hydroxy ketone in which the carbonyl group and the hydroxy group '
                      'are linked by a carbon bearing one hydrogen and one organyl group. Secondary '
                      'alpha-hydroxy ketones are also known as acyloins, and are formally derived '
                      'from reductive coupling of two carboxylic acid groups.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'your_model_name_here',
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
    'success': None,
    'best': None,
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