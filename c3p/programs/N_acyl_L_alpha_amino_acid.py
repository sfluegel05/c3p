"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: CHEBI:58937 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is defined as any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for L-alpha-amino acid backbone: [C@@H](C(=O)O)N
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](C(=O)O)N")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found"

    # Get the nitrogen atom in the amino acid backbone
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    amino_acid_nitrogen = amino_acid_matches[0][2]  # Nitrogen atom index in the amino acid pattern

    # Check for N-acyl substituent: N-C(=O)
    acyl_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No N-acyl substituent found"

    # Get all nitrogen atoms in the molecule
    nitrogen_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if any nitrogen atom is part of an N-acyl group and is connected to the amino acid structure
    for nitrogen in nitrogen_atoms:
        # Check if the nitrogen is part of an N-acyl group
        for neighbor in mol.GetAtomWithIdx(nitrogen).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                # Check if the nitrogen is part of the amino acid structure (backbone or side chain)
                if nitrogen == amino_acid_nitrogen or mol.GetAtomWithIdx(nitrogen).GetDegree() > 1:
                    return True, "Contains L-alpha-amino acid backbone with N-acyl substituent"

    return False, "N-acyl group is not attached to the amino group of the amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58937',
                          'name': 'N-acyl-L-alpha-amino acid',
                          'definition': 'Any L-alpha-amino acid carrying an '
                                        'N-acyl substituent.',
                          'parents': ['CHEBI:58937', 'CHEBI:58937']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48927',
                          'name': 'N-acyl-L-alpha-amino acid',
                          'definition': 'Any L-alpha-amino acid carrying an '
                                        'N-acyl substituent.',
                          'parents': ['CHEBI:21644'],
                          'xrefs': ['KEGG:C02850'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C(CCC)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(N[C@H](C(O)=O)C)=O '
               'NAME: N-arachidonoyl-L-alanine REASON: MISSED N-acyl group is '
               'not attached to the amino group of the amino acid\n'
               ' * SMILES: C(O)([C@@H](N)CCCCNC(C(O)(C)C)=O)=O NAME: '
               'N(6)-(2-hydroxyisobutanoyl)-L-lysine REASON: MISSED N-acyl '
               'group is not attached to the amino group of the amino acid\n'
               ' * SMILES: '
               'C(C(O)=O)C[C@@H](C(=O)O)NC(C1=CC=C(C=C1)S(N2CCC[C@@H]2C(=O)O)(=O)=O)=O '
               'NAME: N-[4-(D-prolinosulfonyl)benzoyl]-L-glutamic acid REASON: '
               'MISSED N-acyl group is not attached to the amino group of the '
               'amino acid\n'
               ' * SMILES: '
               'C=C(C(N[C@@H](CCC([O-])=O)C(=O)[O-])=O)O[C@H]1C=CC=C([C@@H]1O)C([O-])=O '
               'NAME: isochorismoyl-L-glutamate(3-) REASON: MISSED N-acyl '
               'group is not attached to the amino group of the amino acid\n'
               ' * SMILES: CC(=O)N[C@@H](CC(O)=O)C(O)=O NAME: '
               'N-acetyl-L-aspartic acid REASON: MISSED N-acyl group is not '
               'attached to the amino group of the amino acid\n'
               ' * SMILES: '
               'C(C(N[C@H](C(O)=O)CCC(=O)O)=O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: N-[(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl]-L-glutamic '
               'acid REASON: MISSED N-acyl group is not attached to the amino '
               'group of the amino acid\n'
               ' * SMILES: '
               'C1([C@H](CCC1=O)CC(N[C@H](C(O)=O)[C@H](CC)C)=O)C/C=C\\CCO '
               'NAME: N-[(3R)-12-hydroxyjasmonyl]-L-isoleucine REASON: MISSED '
               'N-acyl group is not attached to the amino group of the amino '
               'acid\n'
               ' * SMILES: CC(=O)NCCCC[C@H](N)C(O)=O NAME: '
               'N(6)-acetyl-L-lysine REASON: MISSED N-acyl group is not '
               'attached to the amino group of the amino acid\n'
               ' * SMILES: OC(=O)CC[C@H](NC(=O)CCC(O)=O)C(O)=O NAME: '
               'N(2)-succinyl-L-glutamic acid REASON: MISSED N-acyl group is '
               'not attached to the amino group of the amino acid\n'
               ' * SMILES: NC(=O)N[C@@H](CCC(O)=O)C(O)=O NAME: carglumic acid '
               'REASON: MISSED N-acyl group is not attached to the amino group '
               'of the amino acid\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cl.C=1C=C(OC)C(OC)=C(C1CN2CCNCC2)OC.Cl',
                                     'name': 'Trimetazidine hydrochloride',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3([C@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@]4([H])[C@@H](CCC(NCCS([O-])(=O)=O)=O)C)[H])C)O)[H])C',
                                     'name': 'tauroursocholate',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': '[H]S(=O)(=O)C[C@H](N)C(O)=O',
                                     'name': 'L-cysteine-S-dioxide',
                                     'reason': 'No N-acyl substituent found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C=C(CO[C@@H]1CN(C)CC3=CC=CC=C3OC)N=N2)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'O=C1[C@H]([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@H](O)CC(=C)C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)C',
                                     'name': 'Gilvsin C',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C(=O)C',
                                     'name': '1-[(4S,7S,8S)-5-acetyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-phenylurea',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3=NC=CN=C3',
                                     'name': 'N-[(4S,7R,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[oxo(2-pyrazinyl)methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'COc1ccc(cc1)\\C=C(NC=O)\\C(NC=O)=C\\c1ccc(OS(O)(=O)=O)cc1',
                                     'name': 'Fumiformamide',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'},
                                 {   'smiles': 'CCCCC(CC)COCCCN',
                                     'name': '3-(2-ethylhexoxy)propan-1-amine',
                                     'reason': 'No L-alpha-amino acid backbone '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'The previous code failed because it incorrectly assumed that '
                 'the N-acyl group must be directly attached to the nitrogen '
                 'of the amino acid backbone. However, in many cases, the '
                 'N-acyl group can be attached to other positions on the amino '
                 "acid, such as the side chain (e.g., lysine's epsilon "
                 'nitrogen). The code also did not account for the possibility '
                 'of multiple acyl groups or more complex acyl substituents.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '1. Identify the L-alpha-amino acid backbone.\n'
                 '2. Check for the presence of an N-acyl group (NC(=O)) '
                 'anywhere in the molecule.\n'
                 '3. Ensure that the N-acyl group is attached to a nitrogen '
                 'atom that is part of the amino acid structure (either the '
                 'backbone nitrogen or a side chain nitrogen).\n'
                 '\n'
                 'Here is the revised code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 100,
    'num_true_negatives': 1056,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3103448275862069,
    'recall': 1.0,
    'f1': 0.4736842105263158,
    'accuracy': 0.9167360532889259,
    'negative_predictive_value': 1.0}