"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: CHEBI:36853 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.
    Molecules used for conjugation are: glycine, taurine (and other amino acids); sulfuric acid; glucuronic acid; glucose and other uncharged sugars; and coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for bile acid backbone pattern
    bile_acid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@H]3[C@H]([C@@]4([C@H]([C@@H]5[C@@]4([C@H](C[C@]6([C@@H]([C@@H](O6)C[C@H]2[C@]12C)C)C)C)C)(CC[C@@H]5O)C)C7=CC(=O)C[C@H]7C)[C@H]3C)C"
    if not mol.HasSubstructMatch(bile_acid_pattern):
        return False, "No bile acid backbone found"

    # Look for conjugation patterns
    conjugation_patterns = [
        Chem.MolFromSmarts("[NX3H2,NX4H3+,NX3H0+0]"),  # Amino acids (glycine, taurine, etc.)
        Chem.MolFromSmarts("[SX4](=O)(=O)[O-,OH]"),  # Sulfuric acid/sulfate
        Chem.MolFromSmarts("[OX2H][C@@]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O[CX3](=O)[O-,OH]"),  # Glucuronic acid/glucuronate
        Chem.MolFromSmarts("[OX2H][C@@]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O"),  # Glucose and other uncharged sugars
        Chem.MolFromSmarts("[CX3](=O)[NCCH]") # Coenzyme A
    ]
    conjugation_matches = [mol.GetSubstructMatches(pattern) for pattern in conjugation_patterns]
    if not any(conjugation_matches):
        return False, "No conjugation group found"

    return True, "Contains bile acid backbone with conjugation group"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36853',
        'name': 'bile acid conjugate',
        'definition': 'Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule. Molecules used for conjugation are: glycine, taurine (and other amino acids); sulfuric acid (for which the term ''sulfate'' may be used); glucuronic acid (for which the term ''glucuronate'' may be used); glucose and other uncharged sugars; and coenzyme A.',
        'parents': ['CHEBI:37537', 'CHEBI:36854', 'CHEBI:36962', 'CHEBI:37626', 'CHEBI:37536']
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
    'num_true_positives': 163,
    'num_false_positives': 2,
    'num_true_negatives': 182417,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9878787878787879,
    'recall': 0.9878787878787879,
    'f1': 0.9878787878787879,
    'accuracy': 0.9999906616688736
}