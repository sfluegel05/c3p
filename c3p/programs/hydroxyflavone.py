"""
Classifies: CHEBI:24698 hydroxyflavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxyflavone(smiles: str):
    """
    Determines if a molecule is a hydroxyflavone (flavone with one or more hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains the flavone core structure
    flavone_core_smarts = 'O=C1C=CC(=O)c2c1ccc(O)c2'
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Molecule does not contain flavone core structure"

    # Check for hydroxy groups (-OH) attached to the flavone core
    hydroxy_group_smarts = '[OH]'
    hydroxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts(hydroxy_group_smarts))

    if not hydroxy_groups:
        return False, "No hydroxy groups found"

    return True, "Molecule is a hydroxyflavone with hydroxy groups"

# Example usage
smiles = 'COc1ccc(cc1OC)-c1cc(=O)c2c(O)c(OC)c(O)cc2o1'  # eupatilin
print(is_hydroxyflavone(smiles))  # (True, "Molecule is a hydroxyflavone with hydroxy groups")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24698',
                          'name': 'hydroxyflavone',
                          'definition': 'Any flavone in which one or more ring '
                                        'hydrogens are replaced by hydroxy '
                                        'groups.',
                          'parents': ['CHEBI:24043', 'CHEBI:33822']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 34-35: malformed \\N character escape (<string>, line '
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