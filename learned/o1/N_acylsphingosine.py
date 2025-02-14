"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:51294 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is composed of sphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns

    # Pattern for sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts('[C]-[C@H](O)-[C@H](N)-[C]=[C]-[C]')
    if sphingosine_pattern is None:
        return False, "Invalid sphingosine pattern"

    # Pattern for N-acylation (amide bond)
    n_acyl_pattern = Chem.MolFromSmarts('N-C(=O)-[C]')
    if n_acyl_pattern is None:
        return False, "Invalid N-acylation pattern"

    # Check for sphingosine backbone
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Molecule does not have sphingosine backbone"

    # Check for N-acylation
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "Molecule does not have N-acyl amide bond"

    # Verify the connection between sphingosine backbone and acyl group
    # Find the nitrogen atom in the amide bond
    amide_matches = mol.GetSubstructMatches(n_acyl_pattern)
    sphingosine_matches = mol.GetSubstructMatches(sphingosine_pattern)
    nitrogen_indices_in_amide = [match[0] for match in amide_matches]
    nitrogen_indices_in_sphingosine = [match[2] for match in sphingosine_matches]

    # Check if the nitrogen in the amide bond is part of the sphingosine backbone
    if not set(nitrogen_indices_in_amide).intersection(nitrogen_indices_in_sphingosine):
        return False, "Amide nitrogen is not part of sphingosine backbone"

    return True, "Molecule is an N-acylsphingosine with sphingosine backbone and N-acyl group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51294',
        'name': 'N-acylsphingosine',
        'definition': 'The parent compounds of the ceramide family, composed of sphingosine having an unspecified fatty acyl group attached to the nitrogen.',
        'parents': ['CHEBI:76165', 'CHEBI:76158']
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