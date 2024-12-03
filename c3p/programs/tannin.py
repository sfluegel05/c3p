"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyphenolic structure
    phenol_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetNeighbors()[0].GetSymbol() == 'C':
            if mol.GetBondBetweenAtoms(atom.GetIdx(), atom.GetNeighbors()[0].GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                phenol_count += 1

    if phenol_count < 2:
        return False, "Not enough phenolic hydroxyl groups"

    # Check for glucoside linkage
    glucoside_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 2:
            neighbors = atom.GetNeighbors()
            if neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C':
                if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[0].GetIdx()).GetBondType() == Chem.BondType.SINGLE and \
                   mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[1].GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                    glucoside_found = True
                    break

    if not glucoside_found:
        return False, "No glucoside linkage found"

    return True, "Molecule is a tannin"

# Example usage:
# print(is_tannin('O[C@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4c3c(O)c([C@@H]3[C@@H](O)[C@H](Oc5cc(O)cc(O)c35)c3ccc(O)c(O)c3)c3O[C@@]5(Oc6cc(O)cc(O)c6[C@@H]([C@H]5O)c43)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1'))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26848',
                          'name': 'tannin',
                          'definition': 'Any of a group of astringent '
                                        'polyphenolic vegetable principles or '
                                        'compounds, chiefly complex glucosides '
                                        'of catechol and pyrogallol.',
                          'parents': ['CHEBI:26195']},
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
    'num_true_positives': 21,
    'num_false_positives': 18,
    'num_true_negatives': 2,
    'num_false_negatives': 0,
    'precision': 0.5384615384615384,
    'recall': 1.0,
    'f1': 0.7000000000000001,
    'accuracy': None}