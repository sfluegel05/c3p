"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a quaternary ammonium ion
    quaternary_ammonium_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetFormalCharge() == 1 and atom.GetTotalDegree() == 4:
                quaternary_ammonium_found = True
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:  # Hydrogen atom
                        return False, "Nitrogen has hydrogen atoms attached"
                break

    if quaternary_ammonium_found:
        return True, "Quaternary ammonium ion detected"
    else:
        return False, "No quaternary ammonium ion detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35267',
                          'name': 'quaternary ammonium ion',
                          'definition': 'A derivative of ammonium, NH4(+), in '
                                        'which all four of the hydrogens '
                                        'bonded to nitrogen have been replaced '
                                        'with univalent (usually organyl) '
                                        'groups.',
                          'parents': ['CHEBI:25697', 'CHEBI:35274']},
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
    'num_true_positives': 237,
    'num_false_positives': 14,
    'num_true_negatives': 6,
    'num_false_negatives': 0,
    'precision': 0.9442231075697212,
    'recall': 1.0,
    'f1': 0.9713114754098361,
    'accuracy': None}