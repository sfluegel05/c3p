"""
Classifies: CHEBI:35281 onium betaine
"""
from rdkit import Chem

def is_onium_betaine(smiles: str):
    """
    Determines if a molecule is an onium betaine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an onium betaine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for onium atom (N+, P+, etc.) with no hydrogen atoms
    onium_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0 and atom.GetTotalNumHs() == 0]
    if not onium_atoms:
        return False, "No onium atom with no hydrogen atoms found"

    # Check for anionic atom
    anionic_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0]
    if not anionic_atoms:
        return False, "No anionic atom found"

    # Check that the onium atom is not adjacent to the anionic atom
    for onium_atom in onium_atoms:
        for anionic_atom in anionic_atoms:
            if mol.GetBondBetweenAtoms(onium_atom.GetIdx(), anionic_atom.GetIdx()):
                return False, "Onium atom is adjacent to anionic atom"

    return True, "Molecule is an onium betaine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35281',
                          'name': 'onium betaine',
                          'definition': 'Neutral molecules having '
                                        'charge-separated forms with an onium '
                                        'atom which bears no hydrogen atoms '
                                        'and that is not adjacent to the '
                                        'anionic atom.',
                          'parents': ['CHEBI:27369']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 75-76: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}