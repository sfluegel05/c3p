"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for the presence of an epoxide ring
    epoxide = Chem.MolFromSmarts('C1OC1')
    if not mol.HasSubstructMatch(epoxide):
        return False, "No epoxide ring found"

    # Check for a long carbon chain (fatty acid backbone)
    carbon_chain = Chem.MolFromSmarts('CCCCCCCC')
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chain found"

    return True, "Epoxy fatty acid"

# Example usage
smiles_examples = [
    "CCCCCC\C=C/C\C=C/C\C=C/[C@H](O)[C@H]1O[C@H]1CCCC(O)=O",
    "C(=C\CC1OC1CCCCCO)\CCCCCCCC(=O)O",
    "CCCCCC1OC1C\C=C/C\C=C/C\C=C/CCCC(O)=O",
    "O1[C@H]([C@@H]1[C@H]([C@@H](O)C)C)C[C@@H]2[C@@H](O)[C@@H](O)[C@@H](OC2)C/C(/C)=C/C(O)=O",
    "CCCCC\C=C/C\C=C/C=C/C=C/C[C@H]1O[C@@H]1CCC(O)=O",
    "C(C(O)=O)C/C=C\C/C=C\C/C=C\C/C=C\CC1C(C/C=C\CC)O1",
    "C(CCCO)C/C=C\CC1C(C/C=C\C/C=C\CCCC(O)=O)O1",
    "O=C(CCC/C=C\C/C=C\C/C=C\CC1C(O1)C/C=C\CC)O",
    "CCCC[C@H]1O[C@H]1C\C=C/CCCCCCCC(O)=O",
    "CCCCC\C=C/C[C@H]1O[C@@H]1\C=C\[C@H](O)C\C=C/CCCC(O)=O"
]

for smiles in smiles_examples:
    result, reason = is_epoxy_fatty_acid(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61498',
                          'name': 'epoxy fatty acid',
                          'definition': 'A heterocyclic fatty acid containing '
                                        'an epoxide ring as part of its '
                                        'structure.',
                          'parents': ['CHEBI:23931', 'CHEBI:48847']},
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
              'CCCCCC\\C=C/C\\C=C/C\\C=C/[C@H](O)[C@H]1O[C@H]1CCCC(O)=O -> '
              'True, Reason: Epoxy fatty acid\n'
              'SMILES: C(=C\\CC1OC1CCCCCO)\\CCCCCCCC(=O)O -> True, Reason: '
              'Epoxy fatty acid\n'
              'SMILES: CCCCCC1OC1C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O -> True, '
              'Reason: Epoxy fatty acid\n'
              'SMILES: '
              'O1[C@H]([C@@H]1[C@H]([C@@H](O)C)C)C[C@@H]2[C@@H](O)[C@@H](O)[C@@H](OC2)C/C(/C)=C/C(O)=O '
              '-> True, Reason: Epoxy fatty acid\n'
              'SMILES: CCCCC\\C=C/C\\C=C/C=C/C=C/C[C@H]1O[C@@H]1CCC(O)=O -> '
              'False, Reason: No long carbon chain found\n'
              'SMILES: C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC1C(C/C=C\\CC)O1 '
              '-> False, Reason: No long carbon chain found\n'
              'SMILES: C(CCCO)C/C=C\\CC1C(C/C=C\\C/C=C\\CCCC(O)=O)O1 -> False, '
              'Reason: No long carbon chain found\n'
              'SMILES: O=C(CCC/C=C\\C/C=C\\C/C=C\\CC1C(O1)C/C=C\\CC)O -> '
              'False, Reason: No long carbon chain found\n'
              'SMILES: CCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O -> True, Reason: '
              'Epoxy fatty acid\n'
              'SMILES: '
              'CCCCC\\C=C/C[C@H]1O[C@@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O -> '
              'False, Reason: No long carbon chain found\n',
    'num_true_positives': 5,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 5,
    'precision': 0.7142857142857143,
    'recall': 0.5,
    'f1': 0.588235294117647,
    'accuracy': None}