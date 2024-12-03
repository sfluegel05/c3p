"""
Classifies: CHEBI:35683 aryl sulfide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_aryl_sulfide(smiles: str):
    """
    Determines if a molecule is an aryl sulfide (organic sulfide with sulfur attached to at least one aromatic group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aryl sulfide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfur atoms in the molecule
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if sulfur is attached to at least one aromatic group
    for sulfur in sulfur_atoms:
        neighbors = sulfur.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIsAromatic():
                return True, "Sulfur is attached to at least one aromatic group"

    return False, "Sulfur is not attached to any aromatic group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35683',
                          'name': 'aryl sulfide',
                          'definition': 'Any organic sulfide in which the '
                                        'sulfur is attached to at least one '
                                        'aromatic group.',
                          'parents': ['CHEBI:16385']},
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
    'num_true_positives': 35,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 0,
    'precision': 0.9210526315789473,
    'recall': 1.0,
    'f1': 0.958904109589041,
    'accuracy': None}