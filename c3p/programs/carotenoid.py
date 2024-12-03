"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Carotenoids are tetraterpenoids with C40 skeleton
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 40:
        return False, "Number of carbon atoms is not consistent with C40 skeleton"

    # Check for the presence of conjugated double bonds and other functional groups
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1

    if double_bond_count < 7:
        return False, "Not enough conjugated double bonds"

    # Check for the presence of oxygen-containing functional groups
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if oxygen_count > 4:
        return False, "Too many oxygen-containing functional groups found"

    return True, "Molecule is classified as a carotenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids '
                                        '(C40), formally derived from the '
                                        'acyclic parent, psi,psi-carotene by '
                                        'hydrogenation, dehydrogenation, '
                                        'cyclization, oxidation, or '
                                        'combination of these processes. This '
                                        'class includes carotenes, '
                                        'xanthophylls and certain compounds '
                                        'that arise from rearrangement of the '
                                        'skeleton of psi,psi-carotene or by '
                                        'loss of part of this structure. '
                                        'Retinoids are excluded.',
                          'parents': ['CHEBI:26935']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 10,
    'precision': 1.0,
    'recall': 0.6,
    'f1': 0.7499999999999999,
    'accuracy': None}