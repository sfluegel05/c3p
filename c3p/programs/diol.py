"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol (contains two hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all hydroxy groups (OH)
    hydroxy_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                    hydroxy_groups += 1
                    break

    if hydroxy_groups == 2:
        return True, "Molecule contains two hydroxy groups"
    else:
        return False, f"Molecule contains {hydroxy_groups} hydroxy group(s)"

# Example usage
print(is_diol("OCCCCCCCCCCCCO"))  # 1,12-dodecanediol
print(is_diol("OCC(CO)(C)C"))  # neopentyl glycol
print(is_diol("C1=CC=C(C=C1)O"))  # Phenol (not a diol)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23824',
                          'name': 'diol',
                          'definition': 'A compound that contains two hydroxy '
                                        'groups, generally assumed to be, but '
                                        'not necessarily, alcoholic. Aliphatic '
                                        'diols are also called glycols.',
                          'parents': ['CHEBI:26191']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule contains two hydroxy groups')\n"
              "(True, 'Molecule contains two hydroxy groups')\n"
              "(False, 'Molecule contains 1 hydroxy group(s)')\n",
    'num_true_positives': 29,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 21,
    'precision': 0.9666666666666667,
    'recall': 0.58,
    'f1': 0.725,
    'accuracy': None}