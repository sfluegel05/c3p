"""
Classifies: CHEBI:50584 alkyl alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alkyl_alcohol(smiles: str):
    """
    Determines if a molecule is an alkyl alcohol (an aliphatic alcohol with a hydroxy group at an unspecified position).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkyl alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a hydroxy group (-OH)
    hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            hydroxy = True
            break

    if not hydroxy:
        return False, "No hydroxy group found"

    # Check for the presence of an aliphatic chain (non-aromatic carbon chain)
    aliphatic_chain = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic():
            aliphatic_chain = True
            break

    if not aliphatic_chain:
        return False, "No aliphatic carbon chain found"

    return True, "Molecule is an alkyl alcohol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50584',
                          'name': 'alkyl alcohol',
                          'definition': 'An aliphatic alcohol in which the '
                                        'aliphatic alkane chain is substituted '
                                        'by a hydroxy group at unspecified '
                                        'position.',
                          'parents': ['CHEBI:2571']},
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
    'num_true_positives': 20,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}