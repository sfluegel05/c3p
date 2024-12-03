"""
Classifies: CHEBI:24663 hydroxy-5beta-cholanic acid
"""
from rdkit import Chem

def is_hydroxy_5beta_cholanic_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy-5beta-cholanic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy-5beta-cholanic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one hydroxy group
    hydroxy_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'H' for neigh in atom.GetNeighbors())]
    if not hydroxy_groups:
        return False, "No hydroxy groups found"

    # Check for the presence of the 5beta-cholanic acid core
    # Simplified SMARTS pattern for 5beta-cholanic acid core
    cholanic_acid_core_smarts = "C1CCC2C3CCC4CC[C@@H](C3(C2C1)C)C4"
    core_match = mol.HasSubstructMatch(Chem.MolFromSmarts(cholanic_acid_core_smarts))
    if not core_match:
        return False, "5beta-cholanic acid core not found"

    return True, "Molecule is a hydroxy-5beta-cholanic acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24663',
                          'name': 'hydroxy-5beta-cholanic acid',
                          'definition': 'Any member of the class of '
                                        '5beta-cholanic acids carrying at '
                                        'least one hydroxy group at '
                                        'unspecified position.',
                          'parents': ['CHEBI:33822', 'CHEBI:36248']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 42,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}