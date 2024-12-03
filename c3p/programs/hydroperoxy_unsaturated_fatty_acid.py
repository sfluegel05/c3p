"""
Classifies: CHEBI:194321 hydroperoxy unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_hydroperoxy_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroperoxy unsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroperoxy unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'O' in neighbors and 'O' in neighbors:
                carboxylic_acid = True
                break
    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for at least one C=C double bond
    double_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
            double_bond = True
            break
    if not double_bond:
        return False, "No C=C double bonds found"

    # Check for hydroperoxy group (OOH)
    hydroperoxy = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 2:
                hydroperoxy = True
                break
    if not hydroperoxy:
        return False, "No hydroperoxy group found"

    return True, "Hydroperoxy unsaturated fatty acid detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:194321',
                          'name': 'hydroperoxy unsaturated fatty acid',
                          'definition': 'Any unsaturated fatty acid carrying '
                                        'one or more hydroperoxy substituents.',
                          'parents': ['CHEBI:27208', 'CHEBI:64009']},
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
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}