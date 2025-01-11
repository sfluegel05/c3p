"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: Long-chain fatty acid (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as having a carbon chain length between 13 and 22 carbons,
    with at least one carboxylic acid group (-COOH). The molecule may contain other functional groups
    and cyclic structures as long as the main carbon chain length is within the specified range.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 22:
        return False, f"Carbon chain length {c_count} is not between 13 and 22"

    # Check if the molecule has a reasonable number of rotatable bonds for a fatty acid
    # This helps exclude very rigid structures that are unlikely to be fatty acids
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a long-chain fatty acid"

    # If all checks pass, classify as a long-chain fatty acid
    return True, f"Contains a carboxylic acid group and a carbon chain length of {c_count} (C13 to C22)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15904',
                          'name': 'long-chain fatty acid',
                          'definition': 'A fatty acid with a chain length '
                                        'ranging from C13 to C22.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': ['KEGG:C00638'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: CCCCC\\C=C/CCCC\\C=C/C\\C=C/CCCC(O)=O '
               'NAME: all-cis-icosa-5,8,14-trienoic acid REASON: MISSED Found '
               'unexpected functional group: [OH]\n'
               ' * SMILES: OC(=O)CCCCCCC(CCCCCCCC)C NAME: '
               '8-methyl-hexadecanoic acid REASON: MISSED Found unexpected '
               'functional group: [OH]\n'
               ' * SMILES: C(CCCCCC(=O)O)C/C=C/CC(/C=C(/C=C/CCO)\\O)=O NAME: '
               '(9E,13Z,15E)-14,18-dihydroxy-12-oxooctadeca-9,13,15-trienoic '
               'acid REASON: MISSED Found unexpected functional group: [OH]\n'
               ' * SMILES: OC(=O)CCCCCCCCCCCCC1CCCC1 NAME: '
               'Cyclopentanetridecanoic acid REASON: MISSED Molecule contains '
               'rings, which are not typical for long-chain fatty acids\n'
               ' * SMILES: OCCCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O NAME: '
               '(3R)-3,20-dihydroxyicosanoic acid REASON: MISSED Found '
               'unexpected functional group: [OH]\n'
               ' * SMILES: C(O)(=O)CCCCC/C=C\\CCCCCCCCC(C)C NAME: '
               '17-methyl-7Z-octadecenoic acid REASON: MISSED Found unexpected '
               'functional group: [OH]\n'
               ' * SMILES: OC(=O)CCCCCCC/C=C/C/C=C/CC=C NAME: '
               '9,12,15-hexadecatrienoic acid REASON: MISSED Found unexpected '
               'functional group: [OH]\n'
               ' * SMILES: OC(CCCCC(CC(C)=CC(C)=CC(O)=O)C)C(CO)C(O)=O NAME: '
               '12-Hydroxy-13-(hydroxymethyl)-3,5,7-trimethyltetradeca-2,4-dienedioic '
               'acid REASON: MISSED Multiple carboxylic acid groups found\n'
               ' * SMILES: BrC(C(Br)CCCCCCCC(O)=O)CCCCCCCC NAME: '
               '9,10-Dibromo-stearic acid REASON: MISSED Found unexpected '
               'functional group: [OH]\n'
               ' * SMILES: CCCCC\\C=C/C[C@@H](O)\\C=C\\C=C/C\\C=C/CCCC(O)=O '
               'NAME: 12(R)-HETE REASON: MISSED Found unexpected functional '
               'group: [OH]\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCCC(O)=O',
                                      'name': '10-PAHSA',
                                      'reason': 'Carbon chain length 34 is not '
                                                'between 13 and 22'},
                                  {   'smiles': 'CCCCCCCCCCCCCCC/C(=C(/C(=C(/C(=O)[O-])\\[H])/[H])\\[H])/[H]',
                                      'name': 'Dihomo-linoleate (20:2n6)',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(=O)C/C=C/C=C\\C#CC#CC#CC',
                                      'name': 'Isomycomycin',
                                      'reason': 'Too few rotatable bonds for a '
                                                'long-chain fatty acid'},
                                  {   'smiles': 'C(CCCCC[C@H]1[C@H](C(CC1)=O)C/C=C\\CC)([O-])=O',
                                      'name': '(9R,13R)-1a,1b-dinor-10,11-dihydro-12-oxo-15-phytoenoate',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(C(C(O)/C(=C/C(C(=O)CCCC(CC(=O)N)CC(O)=O)C)/C)C)C(OC)/C=C/CC/C=C/C(O)=O',
                                      'name': '(2E,6E,12E)-19-(2-Amino-2-oxoethyl)-9,11-dihydroxy-8-methoxy-10,12,14-trimethyl-15-oxohenicosa-2,6,12-trienedioic '
                                              'acid',
                                      'reason': 'Carbon chain length 27 is not '
                                                'between 13 and 22'},
                                  {   'smiles': 'O=C(O)C[C@@H](O)/C=C/C=C/C=C/[C@H]([C@@H](O)[C@@H](C[C@H]([C@H](O)CC(=O)C1=CC=C(NC(=O)C)C=C1)C)C)C',
                                      'name': 'Mohangic acid C',
                                      'reason': 'Carbon chain length 28 is not '
                                                'between 13 and 22'},
                                  {   'smiles': 'O=C1[C@@H]([C@@H](CC1)CC(O)=O)C/C=C\\C(O)C',
                                      'name': "Epi-4'-hydroxyjasmonic acid",
                                      'reason': 'Carbon chain length 12 is not '
                                                'between 13 and 22'},
                                  {   'smiles': 'S(O)(=O)(=O)CC(N)(C(OC(=O)CCCCCCCCCCCC(C)C)C(O)/C=C\\CCCCCCCCCCC(C)C)C(O)=O',
                                      'name': "2-carboxy-2-amino-3-O-(13'-methyltetradecanoyl)-4-hydroxy-17-methyloctadec-5-ene-1-sulfonic "
                                              'acid',
                                      'reason': 'Carbon chain length 35 is not '
                                                'between 13 and 22'},
                                  {   'smiles': 'O[C@@](C[N+](C)(C)C)(CC([O-])=O)C=CC=CCCCCCCCCCCO',
                                      'name': 'Hydroxytetradecadienyl-l-carnitine',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O[C@H]1[C@@H]([C@H](C(C1)=C)C/C=C\\CCCC(O)=O)/C=C/[C@@H](O)C(CCCC)(C)C',
                                      'name': '9-Deoxy-9-methylene-16,16-dimethyl '
                                              '-PGE2',
                                      'reason': 'Carbon chain length 23 is not '
                                                'between 13 and 22'}],
    'attempt': 1,
    'reasoning': 'The previous code failed to classify many valid long-chain '
                 'fatty acids because it was too restrictive in its checks. '
                 'Specifically:\n'
                 '\n'
                 '1. It incorrectly rejected molecules with hydroxyl groups '
                 '(-OH), which are common in many natural fatty acids (e.g., '
                 'hydroxy fatty acids).\n'
                 '2. It rejected molecules with rings, even though some fatty '
                 'acids can have cyclic structures (e.g., '
                 'cyclopentane-containing fatty acids).\n'
                 '3. It rejected molecules with multiple carboxylic acid '
                 'groups, which can occur in some fatty acid derivatives.\n'
                 '\n'
                 'To improve the classification, we should:\n'
                 '- Allow hydroxyl groups and other common functional groups '
                 'that can appear in fatty acids.\n'
                 '- Allow cyclic structures, as long as the main carbon chain '
                 'length is within the C13-C22 range.\n'
                 '- Focus primarily on the carbon chain length and the '
                 'presence of at least one carboxylic acid group, while being '
                 'more permissive with other functional groups.\n'
                 '\n'
                 "Here's the revised code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1319,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9295278365045807,
    'f1': 0.9634769905040176,
    'accuracy': 0.9295278365045807,
    'negative_predictive_value': 0.0}