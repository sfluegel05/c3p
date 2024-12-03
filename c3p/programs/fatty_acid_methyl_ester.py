"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic ester group with methanol
    ester_group = Chem.MolFromSmarts('C(=O)OC')
    if not mol.HasSubstructMatch(ester_group):
        return False, "No carboxylic ester group found"

    # Check if the ester group is attached to a long carbon chain (fatty acid)
    ester_matches = mol.GetSubstructMatches(ester_group)
    for match in ester_matches:
        carbonyl_c = match[0]
        ester_o = match[2]
        methanol_c = match[3]

        if mol.GetAtomWithIdx(methanol_c).GetSymbol() == 'C' and mol.GetAtomWithIdx(methanol_c).GetDegree() == 1:
            # Check for a carbon chain attached to the ester group
            carbon_chain = Chem.MolFromSmarts('CCCCCCCC')
            if mol.HasSubstructMatch(carbon_chain):
                return True, "Molecule is a fatty acid methyl ester"
    
    return False, "Molecule does not meet the criteria for fatty acid methyl ester"

# Test cases
test_smiles = [
    "O(C(=O)CCCCCC=CCCCCCCCC)C",
    "O=C(OC)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",
    "O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)C",
    "C(CCCCCCCC)CCCCC(OC)=O",
    "O=C(OC)CC[C@H](CCCCC(C)C)C",
    "O=C1[C@H]([C@@H](CCCC(=O)OC)CC1)C/C=C\\CC",
    "O=C1[C@@H]([C@@H](CCCC(=O)OC)CC1)CCCCC",
    "O(C(=O)CCCCCCCCCCCCCC/C=C/C)C",
    "C(=C\\C(C/C=C\\CCCCCO)OO)/C=C\\C/C=C\\CCCC(=O)OC",
    "O=C(OC)CCCCCC(CCC/C=C/CC)=O",
    "COC(=O)C=C(C)O",
    "CC\\C=C/C\\C=C/C=C/C(O)CCCCCCCC(=O)OC",
    "O=C(OC)CCCCCCCCCC(=CCCCCCC)C",
    "COC(=O)C\\C(=C/C(O)=O)C(O)=O",
    "O1C(C1C/C=C/CC)C/C=C/CCCCCCCC(OC)=O",
    "O=C(O[C@H]1[C@@H]2C=C[C@H]([C@@H](CC)C)[C@@]([C@H]2[C@H](C)C[C@H]1O)(C(=O)CCO)C)/C=C(/CC(=O)OC)\\C",
    "O(C(=O)CC/C=C/CC=C)C"
]

for smiles in test_smiles:
    result, reason = is_fatty_acid_methyl_ester(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:4986',
                          'name': 'fatty acid methyl ester',
                          'definition': 'A fatty acid ester that is the '
                                        'carboxylic ester obtained by the '
                                        'formal condensation of a fatty acid '
                                        'with methanol.',
                          'parents': ['CHEBI:25248', 'CHEBI:35748']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: O(C(=O)CCCCCC=CCCCCCCCC)C, Result: True, Reason: '
              'Molecule is a fatty acid methyl ester\n'
              'SMILES: O=C(OC)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC, Result: '
              'False, Reason: Molecule does not meet the criteria for fatty '
              'acid methyl ester\n'
              'SMILES: O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)C, Result: True, '
              'Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: C(CCCCCCCC)CCCCC(OC)=O, Result: True, Reason: Molecule '
              'is a fatty acid methyl ester\n'
              'SMILES: O=C(OC)CC[C@H](CCCCC(C)C)C, Result: True, Reason: '
              'Molecule is a fatty acid methyl ester\n'
              'SMILES: O=C1[C@H]([C@@H](CCCC(=O)OC)CC1)C/C=C\\CC, Result: '
              'True, Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: O=C1[C@@H]([C@@H](CCCC(=O)OC)CC1)CCCCC, Result: True, '
              'Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: O(C(=O)CCCCCCCCCCCCCC/C=C/C)C, Result: True, Reason: '
              'Molecule is a fatty acid methyl ester\n'
              'SMILES: C(=C\\C(C/C=C\\CCCCCO)OO)/C=C\\C/C=C\\CCCC(=O)OC, '
              'Result: False, Reason: Molecule does not meet the criteria for '
              'fatty acid methyl ester\n'
              'SMILES: O=C(OC)CCCCCC(CCC/C=C/CC)=O, Result: True, Reason: '
              'Molecule is a fatty acid methyl ester\n'
              'SMILES: COC(=O)C=C(C)O, Result: False, Reason: Molecule does '
              'not meet the criteria for fatty acid methyl ester\n'
              'SMILES: CC\\C=C/C\\C=C/C=C/C(O)CCCCCCCC(=O)OC, Result: True, '
              'Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: O=C(OC)CCCCCCCCCC(=CCCCCCC)C, Result: True, Reason: '
              'Molecule is a fatty acid methyl ester\n'
              'SMILES: COC(=O)C\\C(=C/C(O)=O)C(O)=O, Result: False, Reason: '
              'Molecule does not meet the criteria for fatty acid methyl '
              'ester\n'
              'SMILES: O1C(C1C/C=C/CC)C/C=C/CCCCCCCC(OC)=O, Result: True, '
              'Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: '
              'O=C(O[C@H]1[C@@H]2C=C[C@H]([C@@H](CC)C)[C@@]([C@H]2[C@H](C)C[C@H]1O)(C(=O)CCO)C)/C=C(/CC(=O)OC)\\C, '
              'Result: True, Reason: Molecule is a fatty acid methyl ester\n'
              'SMILES: O(C(=O)CC/C=C/CC=C)C, Result: False, Reason: Molecule '
              'does not meet the criteria for fatty acid methyl ester\n',
    'num_true_positives': 12,
    'num_false_positives': 5,
    'num_true_negatives': 12,
    'num_false_negatives': 5,
    'precision': 0.7058823529411765,
    'recall': 0.7058823529411765,
    'f1': 0.7058823529411765,
    'accuracy': None}