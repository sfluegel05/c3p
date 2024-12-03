"""
Classifies: CHEBI:36059 hydroxy monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_hydroxy_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a hydroxy monocarboxylic acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate anion group [C(=O)[O-]]
    carboxylate = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                if any(n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 for n in neighbors) and \
                   any(n.GetAtomicNum() == 8 and n.GetFormalCharge() == 0 and n.GetTotalDegree() == 2 for n in neighbors):
                    carboxylate = True
                    break

    if not carboxylate:
        return False, "No carboxylate anion group found"

    # Check for at least one hydroxy group [-OH]
    hydroxy = any(atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2 and atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())

    if not hydroxy:
        return False, "No hydroxy group found"

    return True, "Molecule is a hydroxy monocarboxylic acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36059',
                          'name': 'hydroxy monocarboxylic acid anion',
                          'definition': 'Any monocarboxylic acid anion '
                                        'carrying at least one hydroxy '
                                        'substituent.',
                          'parents': ['CHEBI:35757']},
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
    'num_false_negatives': 72,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}