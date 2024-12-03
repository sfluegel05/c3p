"""
Classifies: CHEBI:35284 ammonium betaine
"""
from rdkit import Chem

def is_ammonium_betaine(smiles: str):
    """
    Determines if a molecule is an ammonium betaine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ammonium betaine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quaternary ammonium atom (N+ with no hydrogen atoms)
    quaternary_ammonium = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 and atom.GetTotalNumHs() == 0:
            quaternary_ammonium = True
            break

    if not quaternary_ammonium:
        return False, "No quaternary ammonium atom found"

    # Check for adjacent anionic atom
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == -1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() == atom.GetIdx():
                    continue
                if neighbor.GetFormalCharge() == 1:
                    return False, "Anionic atom adjacent to quaternary ammonium atom found"

    return True, "Molecule is an ammonium betaine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35284',
                          'name': 'ammonium betaine',
                          'definition': 'Any neutral molecule having '
                                        'charge-separated forms with a '
                                        'quaternary ammonium atom which bears '
                                        'no hydrogen atoms and that is not '
                                        'adjacent to the anionic atom.',
                          'parents': ['CHEBI:26469', 'CHEBI:35281']},
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
    'num_true_positives': 195,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9948979591836735,
    'recall': 1.0,
    'f1': 0.9974424552429668,
    'accuracy': None}