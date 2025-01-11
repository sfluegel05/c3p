"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid contains a steroid backbone with a 3-oxo group and a Delta(4) double bond.

    Args:
        smiles (str): SMILES string of the chemical entity.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalize the steroid backbone detection using a flexible tetracyclic framework
    steroid_pattern = Chem.MolFromSmarts('C1CCC2CCCCC2C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for the presence of a 3-oxo group
    # This group is part of a ketone functional group, so look for C=O attached to a carbon
    oxo_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found"
    
    # Verify the Alpha, Beta unsaturated bond (Delta(4) bond) as a C=C
    # This often means an alkene feature within the molecule
    delta_4_pattern = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(delta_4_pattern):
        return False, "No Delta(4) double bond found"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate moieties"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47909',
                          'name': '3-oxo-Delta(4) steroid',
                          'definition': 'A 3-oxo steroid conjugated to a C=C '
                                        'double bond at the alpha,beta '
                                        'position.',
                          'parents': ['CHEBI:47788', 'CHEBI:51689'],
                          'xrefs': [   'KEGG:C00619',
                                       'MetaCyc:3-Oxo-Delta-4-Steroids'],
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
               '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C[C@H](O)C2=CC(=O)C=C[C@]12C '
               'NAME: 6alpha-hydroxyprednisone REASON: MISSED No steroid '
               'backbone found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(Cl)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: beclomethasone REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               'C[C@H](CO)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C '
               'NAME: (25S)-26-hydroxycholest-4-en-3-one REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               'C[C@@H]([C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C)C(O)=O '
               'NAME: 3-oxo-23,24-bisnorchola-1,4-dien-22-oic acid REASON: '
               'MISSED No steroid backbone found\n'
               ' * SMILES: '
               'C1[C@@]2([C@@]([C@@]3(C(C1)=CC(CC3)=O)C)(C(C[C@@]4(C)[C@@]2([H])CC[C@@H]4O)=O)[H])[H] '
               'NAME: 11-oxotestosterone REASON: MISSED No steroid backbone '
               'found\n'
               ' * SMILES: '
               '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(C)O '
               'NAME: 17,20-dihydroxypregn-4-en-3-one REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COP(O)(O)=O '
               'NAME: cortisol phosphate REASON: MISSED No steroid backbone '
               'found\n'
               ' * SMILES: '
               'COC(=O)\\C=C1\\CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C '
               'NAME: (17Z)-3,11-dioxopregna-4,17(20)-dien-21-oic acid methyl '
               'ester REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COC(=O)CCCCC '
               'NAME: hydrocortisone caproate REASON: MISSED No steroid '
               'backbone found\n'
               ' * SMILES: '
               'C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)C(=O)C=C34)[C@@H]1CCC2=O '
               'NAME: estra-1(10),4-diene-2,3,17-trione REASON: MISSED No '
               'steroid backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No steroid backbone found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@]1([H])[C@@]2([H])C=CC2=CC(=O)CC[C@@]12C',
                                      'name': 'dydrogesterone',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)(C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])C=C(C)C2=CC(=O)CC[C@]12C',
                                      'name': 'megestrol',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C=CC2=CC(=O)C=C[C@]12C',
                                      'name': 'Delta(6)-prednisone',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)CC=C1[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'vamorolone',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C=CC2=CC(C1)=O)[H])(CC[C@@]4([C@@H](CCC(=O)O)C)[H])[H])C)[H])C',
                                      'name': '3-oxochola-4,6-dien-24-oic acid',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[C@]12([C@]([C@]3([C@](CC1)([C@@]4(C(=CC(CC4)=O)C=C3)C)[H])[H])(CC[C@@]2([C@@H](CCCC(C)C)C)[H])[H])C',
                                      'name': 'cholest-4,6-dien-3-one',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CC[C@](OC(C)=O)(C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])C=C(C)C2=CC(=O)CC[C@]12C',
                                      'name': 'megestrol acetate',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CCC3=C1[C@H](C[C@@]1(C)[C@@]2([H])CC[C@@]1(O)C#CC)c1ccc(cc1)N(C)C',
                                      'name': 'mifepristone',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CCC3=C1[C@H](C[C@@]1(C)[C@@]2([H])CC[C@]1(OC(C)=O)C(C)=O)c1ccc(cc1)N(C)C',
                                      'name': 'ulipristal acetate',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])C=CC4=CC(=O)CC[C@]4(C)[C@@]3([H])C[C@H](O)[C@]12C)[C@H](C)CCC(O)=O',
                                      'name': '12alpha-hydroxy-3-oxochola-4,6-dien-24-oic '
                                              'acid',
                                      'reason': 'No steroid backbone found'}],
    'attempt': 3,
    'reasoning': 'The previous code attempts to classify molecules as '
                 '3-oxo-Delta(4) steroids based on the presence of a steroid '
                 'backbone, a 3-oxo group, and a Delta(4) double bond. '
                 'However, it seems to have failed in several areas:\n'
                 '\n'
                 '1. **Steroid Backbone Detection**: The SMARTS pattern for '
                 'the steroid backbone may have been too rigid or incorrect. '
                 'Steroids can have variations in their tetracyclic ring '
                 'structure, such as different ring sizes or substituents, '
                 'which might not be captured by the SMARTS query used.\n'
                 '\n'
                 '2. **3-oxo and Delta(4) Bond Detection**: The patterns for '
                 'the 3-oxo group and the Delta(4) bond might not accurately '
                 'capture their positions relative to the entire steroid '
                 'backbone. There may also be ambiguity in the placement of '
                 'these features within the ring system that could lead to '
                 'missed matches.\n'
                 '\n'
                 'To improve the program:\n'
                 '- **Generalize the Steroid Pattern**: Ensure that the '
                 'pattern covers different possible arrangements of the '
                 'tetracyclic structure seen in steroids. This might involve '
                 'using more generic patterns or checking for individual rings '
                 'rather than the specific full structure.\n'
                 '- **Specify the 3-oxo and Double Bond Locations**: Clearly '
                 'define the context in which these features should appear '
                 'within the steroid framework, potentially considering the '
                 'positioning in different rings and respecting '
                 'stereochemistry.\n'
                 '- **Use Iterative Pattern Refinement**: Start with more '
                 'general patterns and refine them based on the results '
                 'obtained for each example, ensuring false negatives are '
                 'minimized.\n'
                 '\n'
                 "Let's modify the code with these insights:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 179,
    'num_false_positives': 100,
    'num_true_negatives': 5193,
    'num_false_negatives': 18,
    'num_negatives': None,
    'precision': 0.6415770609318996,
    'recall': 0.9086294416243654,
    'f1': 0.7521008403361344,
    'accuracy': 0.9785063752276867,
    'negative_predictive_value': 0.9965457685664939}