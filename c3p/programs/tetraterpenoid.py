"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C40 skeleton or modified C40 skeleton
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 40:
        return False, f"Number of carbons ({num_carbons}) is less than 40"

    # Check for terpenoid structure
    if not any(atom.GetSymbol() == 'O' for atom in mol.GetAtoms()):
        return False, "No oxygen atoms found, not a terpenoid"

    # Check for multiple double bonds (characteristic of terpenoids)
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 4:  # arbitrary threshold for terpenoid-like structure
        return False, f"Number of double bonds ({num_double_bonds}) is less than 4"

    # Check for typical terpenoid functional groups (e.g., hydroxyl, carbonyl)
    functional_groups = ['O', 'C=O']
    has_functional_group = any(atom.GetSymbol() == 'O' or atom.GetSymbol() == 'C=O' for atom in mol.GetAtoms())
    if not has_functional_group:
        return False, "No typical terpenoid functional groups found"

    return True, "Molecule is a tetraterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26935',
                          'name': 'tetraterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'tetraterpene. The term includes '
                                        'compounds in which the C40 skeleton '
                                        'of the parent tetraterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873']},
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
    'num_true_positives': 22,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 3,
    'precision': 0.9166666666666666,
    'recall': 0.88,
    'f1': 0.8979591836734694,
    'accuracy': None}