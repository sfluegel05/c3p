"""
Classifies: CHEBI:35871 oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors


def is_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is an oxo monocarboxylic acid (Any monocarboxylic acid having at least one additional oxo functional group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for at least one oxo group (C=O) excluding the carboxylic acid group
    oxo_group = Chem.MolFromSmarts('C=O')
    oxo_matches = mol.GetSubstructMatches(oxo_group)
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)

    # Flatten carboxylic matches to get all atoms involved in carboxylic acid groups
    carboxylic_atoms = set(atom for match in carboxylic_matches for atom in match)

    # Check if there is at least one oxo group not part of a carboxylic acid
    oxo_count = 0
    for match in oxo_matches:
        if not any(atom in carboxylic_atoms for atom in match):
            oxo_count += 1

    if oxo_count == 0:
        return False, "No additional oxo group found"

    return True, "Molecule is an oxo monocarboxylic acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35871',
                          'name': 'oxo monocarboxylic acid',
                          'definition': 'Any monocarboxylic acid having at '
                                        'least one additional oxo functional '
                                        'group.',
                          'parents': ['CHEBI:25384', 'CHEBI:25754']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:40:19] SMILES Parse Error: syntax error while parsing: '
             'CC1C\\C(C)=C(C)/C=C/C2=CC=C(C)\\C(O2)=C2OC(C)(C(C(O)=O)C1=O)C(O)=C\x02O\n'
             '[23:40:19] SMILES Parse Error: Failed parsing SMILES '
             "'CC1C\\C(C)=C(C)/C=C/C2=CC=C(C)\\C(O2)=C2OC(C)(C(C(O)=O)C1=O)C(O)=C\x02O' "
             'for input: '
             "'CC1C\\C(C)=C(C)/C=C/C2=CC=C(C)\\C(O2)=C2OC(C)(C(C(O)=O)C1=O)C(O)=C\x02O'\n",
    'stdout': '',
    'num_true_positives': 32,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 2,
    'precision': 0.8,
    'recall': 0.9411764705882353,
    'f1': 0.8648648648648648,
    'accuracy': None}