"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate the molecular weight
    mol_weight = Descriptors.MolWt(mol)
    
    # Check if molecular weight is approximately in the range of sesterterpenoids
    if not (350 <= mol_weight <= 550):
        return False, f"Molecular weight {mol_weight} is not in the typical range for sesterterpenoids (350-550)"

    # Check for the presence of 25 carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 25:
        return False, f"Number of carbon atoms is {carbon_count}, less than 25"

    # Check for terpenoid-like structure (presence of isoprene units)
    isoprene_smarts = Chem.MolFromSmarts("C(C)(C)C(C)C")
    if not mol.HasSubstructMatch(isoprene_smarts):
        return False, "Molecule does not contain isoprene units typically found in terpenoids"

    return True, "Molecule is classified as a sesterterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26660',
                          'name': 'sesterterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'sesterterpene. The term includes '
                                        'compounds in which the C25 skeleton '
                                        'of the parent sesterterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups). Sometimes '
                                        'sesterterpenoids are erroneously '
                                        'referred to as sesterpenoids.',
                          'parents': ['CHEBI:26873']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '[21:05:19] SMILES Parse Error: syntax error while parsing: '
             'O=C/1O[C@H](CC(=O)O)C(\\C1=C\x02/O[C@](CC(CC)C)(C)C=C2)=O\n'
             '[21:05:19] SMILES Parse Error: Failed parsing SMILES '
             "'O=C/1O[C@H](CC(=O)O)C(\\C1=C\x02/O[C@](CC(CC)C)(C)C=C2)=O' for "
             'input: '
             "'O=C/1O[C@H](CC(=O)O)C(\\C1=C\x02/O[C@](CC(CC)C)(C)C=C2)=O'\n",
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.9333333333333333,
    'recall': 0.42424242424242425,
    'f1': 0.5833333333333333,
    'accuracy': None}