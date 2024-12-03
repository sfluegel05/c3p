"""
Classifies: CHEBI:13643 glycol
"""
from rdkit import Chem

def is_glycol(smiles: str):
    """
    Determines if a molecule is a glycol (a diol with hydroxy groups on different carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all hydroxyl groups
    hydroxyls = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'C' for neigh in atom.GetNeighbors())]
    
    if len(hydroxyls) < 2:
        return False, "Less than two hydroxyl groups found"

    # Check if hydroxyl groups are on different carbon atoms
    carbon_atoms_with_oh = [hydroxyl.GetNeighbors()[0] for hydroxyl in hydroxyls if hydroxyl.GetNeighbors()[0].GetSymbol() == 'C']

    if len(carbon_atoms_with_oh) < 2:
        return False, "Hydroxyl groups are not on different carbon atoms"

    # Check if hydroxyl groups are on different carbon atoms
    if len(set(carbon.GetIdx() for carbon in carbon_atoms_with_oh)) < 2:
        return False, "Hydroxyl groups are on the same carbon atom"

    return True, "Valid glycol"

# Example usage
smiles_list = [
    "CCCCCCCCCCCCCCCCCC(O)CCO",
    "COc1cc(ccc1O)C(O)C(CO)c1ccc(O)c(OC)c1",
    "CCCCC\\C=C/C[C@@H](O)C(O)\\C=C\\C(O)C\\C=C/CCCC([O-])=O",
    "OC(CC(O)CCC1=CC(OC)=C(O)C=C1)CCCCC",
    "OC[C@@H](O)COc1ccc(Cl)cc1",
    "OC(CC#CC(O)/C=C/CC)C=C",
    "C[C@H](O)CCO",
    "CCCCCCCCCCCCCCCCCCC[C@H](O)C[C@H](O)CCCC[C@@H](C)[C@H](CC)OC",
    "C=CCCCCCCC(CO)O",
    "C1=2C(=CC(C(C1=O)(OC(C(CC(CC)C)C)=O)C)=O)C=C(OC2)C(C(O)C)O"
]

for smiles in smiles_list:
    result, reason = is_glycol(smiles)
    print(f"SMILES: {smiles}\nIs glycol: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:13643',
                          'name': 'glycol',
                          'definition': 'A diol in which the two hydroxy '
                                        'groups are on different carbon atoms, '
                                        'usually but not necessarily adjacent.',
                          'parents': ['CHEBI:23824']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: CCCCCCCCCCCCCCCCCC(O)CCO\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: COc1cc(ccc1O)C(O)C(CO)c1ccc(O)c(OC)c1\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: '
              'CCCCC\\C=C/C[C@@H](O)C(O)\\C=C\\C(O)C\\C=C/CCCC([O-])=O\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: OC(CC(O)CCC1=CC(OC)=C(O)C=C1)CCCCC\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: OC[C@@H](O)COc1ccc(Cl)cc1\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: OC(CC#CC(O)/C=C/CC)C=C\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: C[C@H](O)CCO\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCC[C@H](O)C[C@H](O)CCCC[C@@H](C)[C@H](CC)OC\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: C=CCCCCCCC(CO)O\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n'
              'SMILES: '
              'C1=2C(=CC(C(C1=O)(OC(C(CC(CC)C)C)=O)C)=O)C=C(OC2)C(C(O)C)O\n'
              'Is glycol: True\n'
              'Reason: Valid glycol\n'
              '\n',
    'num_true_positives': 10,
    'num_false_positives': 10,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}