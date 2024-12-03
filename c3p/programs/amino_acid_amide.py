"""
Classifies: CHEBI:22475 amino acid amide
"""
from rdkit import Chem

def is_amino_acid_amide(smiles: str):
    """
    Determines if a molecule is an amino acid amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxamide group (C(=O)N)
    carboxamide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(carboxamide_pattern):
        return False, "No carboxamide group found"

    # Check for presence of amino acid backbone (NH2-CHR-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Check if carboxyl group (COOH) is converted to carboxamide (C(=O)N)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Carboxyl group not converted to carboxamide"

    return True, "Molecule is an amino acid amide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22475',
                          'name': 'amino acid amide',
                          'definition': 'An amide of an amino acid formed '
                                        'formally by conversion of the carboxy '
                                        'group to a carboxamido group.',
                          'parents': ['CHEBI:37622', 'CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 19-20: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}