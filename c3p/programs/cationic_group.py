"""
Classifies: CHEBI:64766 cationic group
"""
from rdkit import Chem

def is_cationic_group(smiles: str):
    """
    Determines if a molecule belongs to the cationic group (a group that carries an overall positive charge).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cationic group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for overall positive charge
    charge = Chem.GetFormalCharge(mol)
    if charge > 0:
        return True, f"Overall positive charge of {charge}"
    else:
        return False, "No overall positive charge"

# Example usage
print(is_cationic_group("C(=O)([C@@H]([NH3+])CCSC)N[C@H](C(=O)*)[C@H](CC)C"))  # L-methionyl-L-isoleucyl(1+) residue
print(is_cationic_group("CN(C)c1nc2n(c[n+](C)c2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(=O)O-*)[C@@H](O)[C@H]1O"))  # N(2),N(2),N(7)-trimethylguanosine 5'-triphosphate group


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64766',
                          'name': 'cationic group',
                          'definition': 'A group that carries an overall '
                                        'positive charge.',
                          'parents': ['CHEBI:24433']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Overall positive charge of 1')\n"
              "(True, 'Overall positive charge of 1')\n",
    'num_true_positives': 12,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.8571428571428571,
    'f1': 0.923076923076923,
    'accuracy': None}