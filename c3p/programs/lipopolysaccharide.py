"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of trisaccharide repeating unit (two heptose units and octulosonic acid)
    # and 3-hydroxytetradecanoic acid units
    # This is a simplified check based on typical substructures

    # Define substructure for heptose unit (simplified)
    heptose_smarts = '[C@H]1([C@H]([C@H]([C@H]([C@H]([C@H](O1)O)O)O)O)O)O'
    heptose = Chem.MolFromSmarts(heptose_smarts)

    # Define substructure for octulosonic acid (simplified)
    octulosonic_acid_smarts = 'O=C(O)[C@@H]1O[C@H]([C@H]([C@H]([C@H]([C@H]1O)O)O)O)CO'
    octulosonic_acid = Chem.MolFromSmarts(octulosonic_acid_smarts)

    # Define substructure for 3-hydroxytetradecanoic acid (simplified)
    hydroxytetradecanoic_acid_smarts = 'CCCCCCCCCCCCCC(O)C(=O)O'
    hydroxytetradecanoic_acid = Chem.MolFromSmarts(hydroxytetradecanoic_acid_smarts)

    # Check for presence of heptose units
    heptose_matches = mol.GetSubstructMatches(heptose)
    if len(heptose_matches) < 2:
        return False, "Less than two heptose units found"

    # Check for presence of octulosonic acid unit
    octulosonic_acid_matches = mol.GetSubstructMatches(octulosonic_acid)
    if len(octulosonic_acid_matches) < 1:
        return False, "No octulosonic acid unit found"

    # Check for presence of 3-hydroxytetradecanoic acid units
    hydroxytetradecanoic_acid_matches = mol.GetSubstructMatches(hydroxytetradecanoic_acid)
    if len(hydroxytetradecanoic_acid_matches) < 1:
        return False, "No 3-hydroxytetradecanoic acid unit found"

    return True, "Molecule is a lipopolysaccharide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16412',
                          'name': 'lipopolysaccharide',
                          'definition': 'Liposaccharide natural compounds '
                                        'consisting of a trisaccharide '
                                        'repeating unit (two heptose units and '
                                        'octulosonic acid) with '
                                        'oligosaccharide side chains and '
                                        '3-hydroxytetradecanoic acid units '
                                        '(they are a major constituent of the '
                                        'cell walls of Gram-negative '
                                        'bacteria).',
                          'parents': ['CHEBI:35740', 'CHEBI:65212']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}