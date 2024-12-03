"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobases and sugars (ribose and deoxyribose)
    nucleobases = ["C1=NC2=C(N1)N=CN2", "C1=NC2=C(C(=N1)N)N=CN2", "C1=NC2=C(C(=N1)N)N=CO2", "C1=NC=CC(=O)N1", "C1=NC=CC(=O)N1C", "C1=NC=CC(=O)N1C(O)=O"]
    sugars = ["C1C(C(C(O1)CO)O)O", "C1C(C(C(O1)CO)O)O", "C1C(C(C(O1)CO)O)O", "C1C(C(C(O1)CO)O)O", "C1C(C(C(O1)CO)O)O", "C1C(C(C(O1)CO)O)O"]

    # Check for nucleobase
    nucleobase_found = False
    for base in nucleobases:
        base_mol = Chem.MolFromSmarts(base)
        if mol.HasSubstructMatch(base_mol):
            nucleobase_found = True
            break

    if not nucleobase_found:
        return False, "No nucleobase found"

    # Check for sugar (ribose or deoxyribose)
    sugar_found = False
    for sugar in sugars:
        sugar_mol = Chem.MolFromSmarts(sugar)
        if mol.HasSubstructMatch(sugar_mol):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No ribose or deoxyribose found"

    return True, "Molecule is a nucleoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33838',
                          'name': 'nucleoside',
                          'definition': 'An N-glycosyl compound that has both '
                                        'a nucleobase, normally adenine, '
                                        'guanine, xanthine, thymine, cytosine '
                                        'or uracil, and either a ribose or '
                                        'deoxyribose as functional parents.',
                          'parents': [   'CHEBI:21731',
                                         'CHEBI:26912',
                                         'CHEBI:61120']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 57-58: malformed \\N character escape (<string>, line '
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