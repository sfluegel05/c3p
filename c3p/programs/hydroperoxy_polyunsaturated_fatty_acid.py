"""
Classifies: CHEBI:189832 hydroperoxy polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroperoxy_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroperoxy polyunsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroperoxy polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors and neighbors.count('O') == 2:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for at least two double bonds (polyunsaturated)
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds < 2:
        return False, "Less than two double bonds found"

    # Check for hydroperoxy group (OOH)
    hydroperoxy = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 2:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'O' in neighbors:
                hydroperoxy = True
                break

    if not hydroperoxy:
        return False, "No hydroperoxy group found"

    return True, "Molecule is a hydroperoxy polyunsaturated fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:189832',
                          'name': 'hydroperoxy polyunsaturated fatty acid',
                          'definition': 'Any polyunsaturated fatty acid '
                                        'carrying one or more hydroperoxy '
                                        'substituents.',
                          'parents': ['CHEBI:194321', 'CHEBI:26208']},
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
    'num_true_positives': 11,
    'num_false_positives': 1,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.9166666666666666,
    'recall': 1.0,
    'f1': 0.9565217391304348,
    'accuracy': None}