"""
Classifies: CHEBI:37909 oligosaccharide sulfate
"""
from rdkit import Chem

def is_oligosaccharide_sulfate(smiles: str):
    """
    Determines if a molecule is an oligosaccharide sulfate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide sulfate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has at least one O-sulfo substituent (sulfate group)
    sulfate_group = Chem.MolFromSmarts('OS(O)(=O)=O')
    if not mol.HasSubstructMatch(sulfate_group):
        return False, "No sulfate groups found"

    # Check if the molecule is an oligosaccharide
    # Oligosaccharides are made up of 2-10 monosaccharide units
    # Here, we look for the glycosidic bond pattern (C-O-C) to determine the presence of sugar units
    glycosidic_bond = Chem.MolFromSmarts('C-O-C')
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond)

    if len(glycosidic_bond_matches) < 1:
        return False, "No glycosidic bonds found, not an oligosaccharide"

    # If the molecule has both sulfate groups and glycosidic bonds, classify it as an oligosaccharide sulfate
    return True, "Molecule is an oligosaccharide sulfate"

# Example usage:
# smiles = "OS(N[C@H]1[C@H](O[C@@H]([C@H]([C@@H]1OS(O)(=O)=O)O[C@@H]2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O[C@H]3O[C@@H]([C@H]([C@@H]([C@H]3NC(=O)C)O)O[C@@H]4O[C@H]([C@H]([C@@H]([C@H]4OS(O)(=O)=O)O)O[C@@H]5OC(=C[C@@H]([C@H]5OS(O)(=O)=O)O)C(=O)O)COS(O)(=O)=O)C(=O)O)COS(O)(=O)=O)C(=O)O)COS(O)(=O)=O)O)(=O)=O"
# print(is_oligosaccharide_sulfate(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37909',
                          'name': 'oligosaccharide sulfate',
                          'definition': 'Any  carbohydrate sulfate that is an '
                                        'oligosaccharide carrying at least one '
                                        'O-sulfo substituent.',
                          'parents': ['CHEBI:35724', 'CHEBI:63563']},
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
    'num_true_positives': 38,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}