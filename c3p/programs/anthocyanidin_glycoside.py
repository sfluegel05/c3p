"""
Classifies: CHEBI:71583 anthocyanidin glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_glycoside(smiles: str):
    """
    Determines if a molecule is an anthocyanidin glycoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of anthocyanidin core structure
    anthocyanidin_pattern = Chem.MolFromSmarts('c1cc(O)ccc1-c2cc(O)c(O)c(O)c2=[O+]') 
    if not mol.HasSubstructMatch(anthocyanidin_pattern):
        return False, "No anthocyanidin core structure found"

    # Check for glycosyl residues
    glycoside_pattern = Chem.MolFromSmarts('[C@H]1(O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O)')
    matches = mol.GetSubstructMatches(glycoside_pattern)
    if not matches:
        return False, "No glycosyl residues found"

    return True, "Molecule is an anthocyanidin glycoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71583',
                          'name': 'anthocyanidin glycoside',
                          'definition': 'Any anthocyanidin cation having one '
                                        'or more glycosyl residues attached at '
                                        'unspecified positions.',
                          'parents': ['CHEBI:16366', 'CHEBI:24400']},
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
    'num_true_negatives': 16,
    'num_false_negatives': 16,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}