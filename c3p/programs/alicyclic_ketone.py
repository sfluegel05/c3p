"""
Classifies: CHEBI:36132 alicyclic ketone
"""
from rdkit import Chem

def is_alicyclic_ketone(smiles: str):
    """
    Determines if a molecule is an alicyclic ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alicyclic ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a ketone group (C=O)
    ketone = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 3:  # Carbon with 3 bonds (C=O)
                    if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                        ketone = True
                        break
            if ketone:
                break

    if not ketone:
        return False, "No ketone group found"

    # Check if the molecule is cyclic
    if not mol.GetRingInfo().NumRings():
        return False, "Molecule is not cyclic"

    # Check if the molecule is aromatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule is aromatic"

    # Check if the molecule is a benzenoid or contains any aromatic systems
    if any(bond.GetIsAromatic() for bond in mol.GetBonds()):
        return False, "Molecule contains aromatic bonds"

    return True, "Molecule is an alicyclic ketone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36132',
                          'name': 'alicyclic ketone',
                          'definition': 'A cyclic ketone in which the '
                                        'carbocyclic ring structure which may '
                                        'be saturated or unsaturated, but may '
                                        'not be a benzenoid or other aromatic '
                                        'system.',
                          'parents': ['CHEBI:3992']},
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
    'num_true_positives': 75,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 10,
    'precision': 0.9493670886075949,
    'recall': 0.8823529411764706,
    'f1': 0.9146341463414634,
    'accuracy': None}