"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is derived from ammonia where two hydrogen atoms are replaced by hydrocarbyl groups.
    Excludes amide nitrogens and other non-amine nitrogen functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to accurately count hydrogen atoms
    mol = Chem.AddHs(mol)

    # Define improved SMARTS pattern for secondary amine
    # Nitrogen with one hydrogen, degree 3, not bonded to carbonyl carbon, oxygen, nitrogen, or sulfur
    secondary_amine_pattern = Chem.MolFromSmarts("[N;H1;D3;!$(N-C(=O));!$(N-[O,N,S])]")

    # Search for secondary amine substructure
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Molecule contains a secondary amine group"
    else:
        return False, "No secondary amine group found"

__metadata__ = {'chemical_class': {'id': None,
                                   'name': 'secondary amine',
                                   'definition': 'A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.',
                                   'parents': []},
                'config': {'llm_model_name': 'lbl/claude-sonnet',
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
                'accuracy': None}