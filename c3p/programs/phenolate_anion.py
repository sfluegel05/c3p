"""
Classifies: CHEBI:50525 phenolate anion
"""
from rdkit import Chem

def is_phenolate_anion(smiles: str):
    """
    Determines if a molecule is a phenolate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenolate anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the phenolate anion group
    phenolate_anion_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == -1:
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                    phenolate_anion_found = True
                    break
            if phenolate_anion_found:
                break

    if not phenolate_anion_found:
        return False, "No phenolate anion group found"

    return True, "Phenolate anion group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50525',
                          'name': 'phenolate anion',
                          'definition': 'An organic anion arising from '
                                        'deprotonation of the OH function of a '
                                        'phenol compound.',
                          'parents': ['CHEBI:25696']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 141-142: malformed \\N character escape (<string>, line '
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