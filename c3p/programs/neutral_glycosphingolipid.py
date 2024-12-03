"""
Classifies: CHEBI:25513 neutral glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neutral_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a neutral glycosphingolipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neutral glycosphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosphingolipid core structure
    # Glycosphingolipids typically have a long-chain base (sphingosine or similar) and one or more sugar residues
    found_sphingosine = False
    found_sugar = False

    # Identify sphingosine-like structure
    sphingosine_smarts = Chem.MolFromSmarts("C[C@@H](O)[C@H](NC=O)CO")
    if mol.HasSubstructMatch(sphingosine_smarts):
        found_sphingosine = True

    # Identify sugar residues (simple heuristic: presence of multiple hydroxyl groups)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            hydroxyl_count += 1
    if hydroxyl_count > 3:
        found_sugar = True

    if not found_sphingosine:
        return False, "No sphingosine-like structure found"
    if not found_sugar:
        return False, "No sugar residues found"

    # Check if glycosyl moieties are unsubstituted
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'S' for neighbor in atom.GetNeighbors()):
            return False, "Glycosyl moieties are substituted"

    return True, "Neutral glycosphingolipid"

# Test examples
smiles_examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\C=C\CCCCCCCCCCCCC",
    "CCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\C=C\CCCCCCCCCCCCC",
    "CCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H]2O)[C@H](O)[C@H]1O)NC([*])=O",
    "O[C@@H]1[C@H]([C@H](O[C@@H]([C@@H]1O)CO)O[C@H]2[C@@H]([C@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@@H]([C@H](O[C@@H]([C@@H]3O)CO)O[C@H]4[C@@H]([C@H](O[C@@H]([C@@H]4O)CO)O[C@H]5[C@@H]([C@H](O[C@@H]([C@@H]5O)CO)O[C@H]6[C@@H]([C@H](O[C@@H]([C@@H]6O)CO)OC[C@@H]([C@@H](*)O)NC(=O)*)O)O)O)NC(=O)C)O)O"
]

for smiles in smiles_examples:
    result, reason = is_neutral_glycosphingolipid(smiles)
    print(f"SMILES: {smiles}\nIs Neutral Glycosphingolipid: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25513',
                          'name': 'neutral glycosphingolipid',
                          'definition': 'Any glycosphingolipid containing '
                                        'unsubstituted glycosyl moieties.',
                          'parents': ['CHEBI:24402']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC\n'
              'Is Neutral Glycosphingolipid: True\n'
              'Reason: Neutral glycosphingolipid\n'
              '\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC\n'
              'Is Neutral Glycosphingolipid: True\n'
              'Reason: Neutral glycosphingolipid\n'
              '\n'
              'SMILES: '
              'CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H]2O)[C@H](O)[C@H]1O)NC([*])=O\n'
              'Is Neutral Glycosphingolipid: True\n'
              'Reason: Neutral glycosphingolipid\n'
              '\n'
              'SMILES: '
              'O[C@@H]1[C@H]([C@H](O[C@@H]([C@@H]1O)CO)O[C@H]2[C@@H]([C@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@@H]([C@H](O[C@@H]([C@@H]3O)CO)O[C@H]4[C@@H]([C@H](O[C@@H]([C@@H]4O)CO)O[C@H]5[C@@H]([C@H](O[C@@H]([C@@H]5O)CO)O[C@H]6[C@@H]([C@H](O[C@@H]([C@@H]6O)CO)OC[C@@H]([C@@H](*)O)NC(=O)*)O)O)O)NC(=O)C)O)O\n'
              'Is Neutral Glycosphingolipid: True\n'
              'Reason: Neutral glycosphingolipid\n'
              '\n',
    'num_true_positives': 64,
    'num_false_positives': 13,
    'num_true_negatives': 7,
    'num_false_negatives': 11,
    'precision': 0.8311688311688312,
    'recall': 0.8533333333333334,
    'f1': 0.8421052631578949,
    'accuracy': None}