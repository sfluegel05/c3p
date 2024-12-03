"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound (contains at least one carbon-bromine bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Check if the bond is between a carbon and a bromine atom
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'Br') or (begin_atom.GetSymbol() == 'Br' and end_atom.GetSymbol() == 'C'):
            return True, "Contains at least one carbon-bromine bond"
    
    return False, "No carbon-bromine bonds found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37141',
                          'name': 'organobromine compound',
                          'definition': 'A compound containing at least one '
                                        'carbon-bromine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:22928']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 12-13: malformed \\N character escape (<string>, line '
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