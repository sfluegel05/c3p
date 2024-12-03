"""
Classifies: CHEBI:79323 4-hydroxybenzoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_4_hydroxybenzoate_ester(smiles: str):
    """
    Determines if a molecule is a 4-hydroxybenzoate ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4-hydroxybenzoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 4-hydroxybenzoic acid core substructure
    hydroxybenzoic_acid_core = Chem.MolFromSmarts('Oc1ccc(C(=O)O)cc1')
    if not mol.HasSubstructMatch(hydroxybenzoic_acid_core):
        return False, "No 4-hydroxybenzoic acid core found"

    # Define the ester linkage substructure
    ester_linkage = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(ester_linkage):
        return False, "No ester linkage found"

    return True, "4-hydroxybenzoate ester found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:79323',
                          'name': '4-hydroxybenzoate ester',
                          'definition': 'A benzoate ester that is an ester of '
                                        '4-hydroxybenzoic acid.',
                          'parents': ['CHEBI:33853', 'CHEBI:36054']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}