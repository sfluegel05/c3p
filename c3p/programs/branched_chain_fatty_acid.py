"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of carboxylic acid group
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)O')):
        return False, "No carboxylic acid group found"

    # Check for presence of alkyl substituents on the hydrocarbon chain
    chain_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetDegree() == 2]
    branches = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetDegree() > 2]

    if not branches:
        return False, "No alkyl substituents found"
    
    # Check if the alkyl substituents are typical (methyl, ethyl, etc.)
    for branch in branches:
        if any(neighbor.GetSymbol() not in ['C', 'H', 'O'] for neighbor in branch.GetNeighbors()):
            return False, "Branch contains non-carbon atoms"
    
    return True, "Molecule is a branched-chain fatty acid"

# Example usage
smiles_examples = [
    "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",
    "CCCCCCCCCCCCCCCCCCCCCCC[C@H]([C@H](O)CCCCCCCCCCCCCCCCC[C@H]1C[C@H]1CCCCCCCCCCCCCCCC[C@H](OC)[C@@H](C)CCCCCCCCCCCCCCCCCC)C(O)=O",
    "OC(=O)C(CCC)/C=C/C",
    "CC(C)CCCCCCCCCC(O)=O",
    "OC(=O)C(CCC)(CC)C",
    "C(CCCCCC1C(C1)CCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)CCCCCCCCC2C(C2)CCCCCCCCCCCCCCCCCCCC",
    "C(=O)(O)*",
    "C(=O)(O)*",
    "[H]\\C(C)=C(\\C)C(O)=O",
    "OC(=O)C/C(/C)=C\\C",
    "[H]C(C)=C(C)C(O)=O",
    "CCCCCCCCCCCCCC(C)C(O)=O",
    "OC(=O)C(CCC)CC",
    "OC(=O)C(C(NC(OC)=O)(C)C)(C)C",
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)C(CCCCCCCCCCCCCC)C(O)=O",
    "CCCC(C)C(O)=O",
    "CCC(C)CCCCCCCCC(O)=O",
    "CCC(C)CCCCCCCCCCC(O)=O",
    "C(CCCCCC(CCC)C)CCCCC(O)=O",
    "CCCCCCCCCCCCCCCC(O)C(CCCCCCCCCCCCCC)C(O)=O",
    "O=C(O)[C@H](C(=C)C)CCC(=O)OC",
    "C(CCCCCC1C(C1)CCCCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)CCCCCCCC/C=C\\CCCCCCCCCCCCCCCCCCCC"
]

for smiles in smiles_examples:
    result, reason = is_branched_chain_fatty_acid(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35819',
                          'name': 'branched-chain fatty acid',
                          'definition': 'Any fatty acid in which the parent '
                                        'hydrocarbon chain has one or more '
                                        'alkyl substituents; a common '
                                        'component in animal and bacterial '
                                        'lipids. The fatty acyl chain is '
                                        'usually saturated and the substituent '
                                        'a methyl group; however, unsaturated '
                                        'BCFAs are found in marine animals, '
                                        'and branches other than methyl are '
                                        'found in microbial lipids.',
                          'parents': ['CHEBI:35366']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O -> True, Reason: '
              'Molecule is a branched-chain fatty acid\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCC[C@H]([C@H](O)CCCCCCCCCCCCCCCCC[C@H]1C[C@H]1CCCCCCCCCCCCCCCC[C@H](OC)[C@@H](C)CCCCCCCCCCCCCCCCCC)C(O)=O '
              '-> True, Reason: Molecule is a branched-chain fatty acid\n'
              'SMILES: OC(=O)C(CCC)/C=C/C -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: CC(C)CCCCCCCCCC(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: OC(=O)C(CCC)(CC)C -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: '
              'C(CCCCCC1C(C1)CCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)CCCCCCCCC2C(C2)CCCCCCCCCCCCCCCCCCCC '
              '-> True, Reason: Molecule is a branched-chain fatty acid\n'
              'SMILES: C(=O)(O)* -> False, Reason: Branch contains non-carbon '
              'atoms\n'
              'SMILES: C(=O)(O)* -> False, Reason: Branch contains non-carbon '
              'atoms\n'
              'SMILES: [H]\\C(C)=C(\\C)C(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: OC(=O)C/C(/C)=C\\C -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: [H]C(C)=C(C)C(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: CCCCCCCCCCCCCC(C)C(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: OC(=O)C(CCC)CC -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: OC(=O)C(C(NC(OC)=O)(C)C)(C)C -> False, Reason: Branch '
              'contains non-carbon atoms\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)C(CCCCCCCCCCCCCC)C(O)=O -> '
              'True, Reason: Molecule is a branched-chain fatty acid\n'
              'SMILES: CCCC(C)C(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: CCC(C)CCCCCCCCC(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: CCC(C)CCCCCCCCCCC(O)=O -> True, Reason: Molecule is a '
              'branched-chain fatty acid\n'
              'SMILES: C(CCCCCC(CCC)C)CCCCC(O)=O -> True, Reason: Molecule is '
              'a branched-chain fatty acid\n'
              'SMILES: CCCCCCCCCCCCCCCC(O)C(CCCCCCCCCCCCCC)C(O)=O -> True, '
              'Reason: Molecule is a branched-chain fatty acid\n'
              'SMILES: O=C(O)[C@H](C(=C)C)CCC(=O)OC -> True, Reason: Molecule '
              'is a branched-chain fatty acid\n'
              'SMILES: '
              'C(CCCCCC1C(C1)CCCCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)CCCCCCCC/C=C\\CCCCCCCCCCCCCCCCCCCC '
              '-> True, Reason: Molecule is a branched-chain fatty acid\n',
    'num_true_positives': 19,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 3,
    'precision': 0.48717948717948717,
    'recall': 0.8636363636363636,
    'f1': 0.6229508196721312,
    'accuracy': None}