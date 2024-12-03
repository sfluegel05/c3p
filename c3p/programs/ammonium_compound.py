"""
Classifies: CHEBI:35276 ammonium compound
"""
from rdkit import Chem

def is_ammonium_compound(smiles: str):
    """
    Determines if a molecule is an ammonium compound (NH4(+))Y(-) and derivatives, in which one or more of the hydrogens bonded to nitrogen have been replaced with univalent groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ammonium compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ammonium ion [NH4+]
    ammonium_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
            hydrogen_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'H')
            if hydrogen_count >= 1:
                ammonium_found = True
                break

    if not ammonium_found:
        return False, "No ammonium ion found"

    # Check for univalent group replacements on nitrogen
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
            hydrogen_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'H')
            if hydrogen_count < 4:
                return True, "Ammonium compound with univalent group replacements found"
    
    return False, "No univalent group replacements found on ammonium ion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35276',
                          'name': 'ammonium compound',
                          'definition': 'Compounds (NH4(+))Y(-) and '
                                        'derivatives, in which one or more of '
                                        'the hydrogens bonded to nitrogen have '
                                        'been replaced with univalent groups.',
                          'parents': ['CHEBI:51143']},
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
    'num_false_negatives': 25,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}