"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amide group
    amide_group = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_group):
        return False, "No amide group found"

    # Check for the presence of a sphingoid base
    sphingoid_base = Chem.MolFromSmarts("C[C@H](O)C[NH]")
    if not mol.HasSubstructMatch(sphingoid_base):
        return False, "No sphingoid base found"

    # Check for the presence of a long-chain fatty acid (14-26 carbons)
    fatty_acid = Chem.MolFromSmarts("C(=O)C" + "C" * 13)
    if not mol.HasSubstructMatch(fatty_acid):
        return False, "No long-chain fatty acid found"

    # Check for the presence of hydroxyl group on carbon 2
    hydroxyl_group = Chem.MolFromSmarts("C[C@H](O)C")
    if not mol.HasSubstructMatch(hydroxyl_group):
        return False, "No hydroxyl group on carbon 2 found"

    return True, "Molecule is a ceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17761',
                          'name': 'ceramide',
                          'definition': 'Ceramides (N-acyl-sphingoid bases) '
                                        'are a major subclass of sphingoid '
                                        'base derivatives with an amide-linked '
                                        'fatty acid. The fatty acids are '
                                        'typically saturated or '
                                        'monounsaturated with chain lengths '
                                        'from 14 to 26 carbon atoms; the '
                                        'presence of a hydroxyl group on '
                                        'carbon 2 is fairly common. Ceramides '
                                        'are generally precursors of more '
                                        'complex sphingolipids. In the '
                                        'illustrated generalised structure, '
                                        'R(1) = OH, OX (where X = acyl, '
                                        'glycosyl, phosphate, phosphonate, '
                                        'etc.), or H.',
                          'parents': ['CHEBI:26739', 'CHEBI:37622']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 59,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 101,
    'precision': 0.9365079365079365,
    'recall': 0.36875,
    'f1': 0.5291479820627802,
    'accuracy': None}