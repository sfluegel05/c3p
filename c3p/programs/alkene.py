"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Check if the molecule is acyclic
    if not rdMolDescriptors.CalcNumRings(mol) == 0:
        return False, "Molecule contains rings"

    # Count the number of carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and \
           bond.GetBeginAtom().GetAtomicNum() == 6 and \
           bond.GetEndAtom().GetAtomicNum() == 6:
            double_bond_count += 1

    if double_bond_count != 1:
        return False, f"Found {double_bond_count} carbon-carbon double bonds, need exactly 1"

    # Check if the molecule is a hydrocarbon (only contains C and H)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            return False, "Molecule contains non-carbon and non-hydrogen atoms"

    # Check the general formula CnH2n
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if h_count != 2 * c_count:
        return False, f"Hydrogen count ({h_count}) does not match 2 * carbon count ({2 * c_count})"

    return True, "Acyclic hydrocarbon with exactly one carbon-carbon double bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32878',
                          'name': 'alkene',
                          'definition': 'An acyclic branched or unbranched '
                                        'hydrocarbon having one carbon-carbon '
                                        'double bond and the general formula '
                                        'CnH2n. Acyclic branched or unbranched '
                                        'hydrocarbons having more than one '
                                        'double bond are alkadienes, '
                                        'alkatrienes, etc.',
                          'parents': ['CHEBI:32877', 'CHEBI:33281']},
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


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32878',
                          'name': 'alkene',
                          'definition': 'An acyclic branched or unbranched '
                                        'hydrocarbon having one carbon-carbon '
                                        'double bond and the general formula '
                                        'CnH2n. Acyclic branched or unbranched '
                                        'hydrocarbons having more than one '
                                        'double bond are alkadienes, '
                                        'alkatrienes, etc.',
                          'parents': ['CHEBI:33645'],
                          'xrefs': ['KEGG:C01372'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: CCCCCCCCC=C NAME: 1-decene REASON: '
               'MISSED Hydrogen count (0) does not match 2 * carbon count '
               '(20)\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCC=C NAME: nonadec-1-ene REASON: '
               'MISSED Hydrogen count (0) does not match 2 * carbon count '
               '(38)\n'
               ' * SMILES: CCCCCCCCCCCCCCC=C NAME: 1-hexadecene REASON: MISSED '
               'Hydrogen count (0) does not match 2 * carbon count (32)\n'
               ' * SMILES: CCCCC\\C=C/C NAME: (Z)-2-octene REASON: MISSED '
               'Hydrogen count (0) does not match 2 * carbon count (16)\n'
               ' * SMILES: CCCCCC(C)CCC=C NAME: 5-methyldec-1-ene REASON: '
               'MISSED Hydrogen count (0) does not match 2 * carbon count '
               '(22)\n'
               ' * SMILES: CCCCCCCCCCC=C NAME: 1-dodecene REASON: MISSED '
               'Hydrogen count (0) does not match 2 * carbon count (24)\n'
               ' * SMILES: C=C NAME: ethene REASON: MISSED Hydrogen count (0) '
               'does not match 2 * carbon count (4)\n'
               ' * SMILES: [H]C(C)=C([H])CCCCC NAME: 2-octene REASON: MISSED '
               'Hydrogen count (0) does not match 2 * carbon count (16)\n'
               ' * SMILES: CC\\C=C/CC NAME: cis-3-hexene REASON: MISSED '
               'Hydrogen count (0) does not match 2 * carbon count (12)\n'
               ' * SMILES: CCCCCCCCCC\\C=C\\C NAME: (2E)-2-tridecene REASON: '
               'MISSED Hydrogen count (0) does not match 2 * carbon count '
               '(26)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'Found 0 carbon-carbon double '
                                               'bonds, need exactly 1'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'Molecule contains rings'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'Molecule contains rings'}],
    'sample_false_negatives': [   {   'smiles': 'C(CCCC(CCCC(C=C)=C)C)(CCCC(C)C)C',
                                      'name': 'neophytadiene',
                                      'reason': 'Found 2 carbon-carbon double '
                                                'bonds, need exactly 1'},
                                  {   'smiles': 'C=CCCCCCCC1OC1',
                                      'name': '1,2-epoxy-9-decene',
                                      'reason': 'Molecule contains rings'}],
    'attempt': 1,
    'reasoning': 'The previous program failed because it incorrectly '
                 'calculated the hydrogen count in the SMILES strings. The '
                 "issue arises because RDKit's `MolFromSmiles` function does "
                 'not explicitly include hydrogens in the molecular graph '
                 'unless they are explicitly specified in the SMILES string. '
                 'This leads to an incorrect hydrogen count, which in turn '
                 'causes the program to fail the `CnH2n` check.\n'
                 '\n'
                 'To fix this, we need to add hydrogens to the molecule '
                 'explicitly before counting the atoms. This can be done using '
                 'the `AddHs` function from RDKit. Additionally, we should '
                 'ensure that the molecule is acyclic and contains exactly one '
                 'carbon-carbon double bond, as before.\n'
                 '\n'
                 'Here is the corrected code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 57,
    'num_false_positives': 39,
    'num_true_negatives': 142202,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.59375,
    'recall': 0.9661016949152542,
    'f1': 0.7354838709677419,
    'accuracy': 0.9997118763176388,
    'negative_predictive_value': 0.9999859356980113}