"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.

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

    # Define amine SMARTS pattern
    amine_pattern = Chem.MolFromSmarts('[NX3;!$(N-C=O);!$(N-C#N);!$(N=O);!$(N=N);!$(N=*);!$([N+])]')

    # Search for amine groups
    if mol.HasSubstructMatch(amine_pattern):
        num_amines = len(mol.GetSubstructMatches(amine_pattern))
        return True, f"Found {num_amines} amine group(s)"
    else:
        return False, "No amine groups found"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32952',
                                         'name': 'amine',
                                         'definition': 'A compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.',
                                         'parents': [] },
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