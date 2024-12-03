"""
Classifies: CHEBI:131871 hydroxy polyunsaturated fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_polyunsaturated_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a hydroxy polyunsaturated fatty acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy polyunsaturated fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for an anion group (carboxylate group)
    carboxylate_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == -1:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C':
                carboxylate_group = True
                break

    if not carboxylate_group:
        return False, "No carboxylate group found"

    # Check for polyunsaturation (multiple double bonds)
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds < 2:
        return False, "Not enough double bonds for polyunsaturation"

    # Check for hydroxy groups
    hydroxy_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == 0:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C':
                hydroxy_groups += 1

    if hydroxy_groups < 1:
        return False, "No hydroxy groups found"

    return True, "Molecule is a hydroxy polyunsaturated fatty acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:131871',
                          'name': 'hydroxy polyunsaturated fatty acid anion',
                          'definition': 'Any polyunsaturated fatty acid anion '
                                        'carrying one or more hydroxy '
                                        'substituents.',
                          'parents': ['CHEBI:59835', 'CHEBI:76567']},
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
    'num_true_positives': 12,
    'num_false_positives': 5,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7058823529411765,
    'recall': 1.0,
    'f1': 0.8275862068965517,
    'accuracy': None}