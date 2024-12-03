"""
Classifies: CHEBI:71939 diterpene glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_diterpene_glycoside(smiles: str):
    """
    Determines if a molecule is a diterpene glycoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpene glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has a glycoside moiety
    glycoside_substructure = Chem.MolFromSmarts('[C@H]1(O[C@@H]([C@H](O)[C@@H](O)[C@H]1O)CO)')
    if not mol.HasSubstructMatch(glycoside_substructure):
        return False, "No glycoside moiety found"

    # Check if molecule has a diterpenoid moiety
    diterpenoid_substructure = Chem.MolFromSmarts('C1(CC[C@@H](C(C1)C)C)C')
    if not mol.HasSubstructMatch(diterpenoid_substructure):
        return False, "No diterpenoid moiety found"

    return True, "Molecule is a diterpene glycoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71939',
                          'name': 'diterpene glycoside',
                          'definition': 'A terpene glycoside in which the '
                                        'terpene moiety is a diterpenoid.',
                          'parents': ['CHEBI:23849', 'CHEBI:61777']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 37-38: malformed \\N character escape (<string>, line '
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