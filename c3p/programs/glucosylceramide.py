"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide has a sphingosine backbone with an amide-linked fatty acyl chain
    and a β-D-glucose sugar attached to the primary hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for β-D-glucose moiety
    # This pattern explicitly matches β-D-glucosyl linkage
    glucose_pattern = Chem.MolFromSmarts('C1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    
    # Check for the glucosyl moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found"
    
    # Define a sphingosine-like backbone pattern with amine and hydroxyl groups
    sphingosine_backbone_pattern = Chem.MolFromSmarts('N[C@@H](CO)C')
    
    # Ensure that the substructure for a sphingosine-like backbone is present
    if not mol.HasSubstructMatch(sphingosine_backbone_pattern):
        return False, "No recognizable sphingosine backbone pattern found"
    
    # Amide linkage pattern for fatty acyl chain attachment
    amide_pattern = Chem.MolFromSmarts('C(=O)N[C@@H]')
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Amide linkage characteristic of fatty acyl attachment not found"
    
    # If all essential patterns match, classify as glucosylceramide
    return True, "Molecule contains a sphingosine backbone, amide linkage, and β-D-glucose moiety"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36500',
                          'name': 'glucosylceramide',
                          'definition': 'Any of the cerebrosides in which the '
                                        'monosaccharide head group is glucose.',
                          'parents': ['CHEBI:23079', 'CHEBI:62941'],
                          'xrefs': ['PMID:16758576'],
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
               'C(CCCCCC(CC)C)CC\\C=C\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCCCCCCCCC)=O)CO[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO '
               'NAME: '
               "beta-D-glucosyl-(1<->1')-N-pentacosanoyl-14-methylhexadecasphingosine "
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               '[C@H]([C@@H](/C=C/CCCCCCCC/C=C\\CCC)O)(NC(CCCCCCCCCCCCCCCCCCCC)=O)CO[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)CO '
               'NAME: '
               "beta-D-glucosyl-(1<->1')-N-henicosanoyl-(4E,14Z)-sphingadienine "
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxyheptadecanoyl)-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'C(CCCCCC(CC)C)CC\\C=C\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCCCCCCC)=O)CO[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO '
               'NAME: '
               "beta-D-glucosyl-(1<->1')-N-tricosanoyl-14-methylhexadecasphingosine "
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: beta-D-glucosyl-N-octadecanoylsphingosine REASON: MISSED '
               'No β-D-glucose moiety found attached to the primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxytetracosanoyl)-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: beta-D-glucosyl-N-(hexacosanoyl)sphingosine REASON: '
               'MISSED No β-D-glucose moiety found attached to the primary '
               'hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCC(C)C '
               'NAME: '
               'N-octacosanoyl-1-O-beta-D-glucosyl-15-methylhexadecasphing-4-enine '
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCC(C)C '
               'NAME: '
               'N-(2-hydroxyhenicosanoyl)-1-O-beta-D-glucosyl-4-hydroxy-15-methylhexadecasphinganine '
               'REASON: MISSED No β-D-glucose moiety found attached to the '
               'primary hydroxyl\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCC\\C=C\\[C@@H](O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CC\\C=C(/C)CCCCCC '
               'NAME: chrysogeside B REASON: MISSED No β-D-glucose moiety '
               'found attached to the primary hydroxyl\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No recognizable sphingosine '
                                               'backbone pattern found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No recognizable sphingosine '
                                               'backbone pattern found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No β-D-glucose moiety found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No β-D-glucose moiety found'}],
    'sample_false_negatives': [   {   'smiles': '[C@H]1([C@H](O[C@@H](OC[C@@H]([C@@H](CCCCCCCCCCCCCCC)O)OC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@@H]([C@H]1O)O)CO)O',
                                      'name': "beta-D-glucosyl-(1<->1')-N-[(15Z)-tetracosenoyl]sphinganine",
                                      'reason': 'No recognizable sphingosine '
                                                'backbone pattern found'}],
    'attempt': 4,
    'reasoning': 'To classify a molecule as a glucosylceramide, the essential '
                 'components are:\n'
                 '\n'
                 '1. A sphingosine backbone, which includes a long-chain base '
                 'with an amine group.\n'
                 '2. An amide-linked fatty acyl chain.\n'
                 '3. A β-D-glucose moiety attached to the primary hydroxyl '
                 'group.\n'
                 '\n'
                 'The previous implementation failed primarily due to the '
                 'incorrect detection of the β-D-glucose moiety, which '
                 'resulted in all false negatives as reported. This failure '
                 'likely stems from the incorrect or too stringent definition '
                 'of the β-D-glucose substructure pattern. Here are some '
                 'specific improvements:\n'
                 '\n'
                 '- **Sphingosine backbone detection**: Ensure that the '
                 'sphingosine pattern captures typical arrangements, including '
                 'the hydroxyl and amine groups.\n'
                 '- **β-D-glucose moiety**: The pattern must accurately '
                 'reflect the stereochemistry and linkage to the sphingosine.\n'
                 '- **General structure**: Double-check the linkage patterns '
                 'and simplify or refine where necessary without being too '
                 'broad.\n'
                 '\n'
                 'We can revise the program to include a more comprehensive '
                 'set of patterns and explanations for detected features:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 79,
    'num_false_positives': 100,
    'num_true_negatives': 1292,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.441340782122905,
    'recall': 0.9875,
    'f1': 0.61003861003861,
    'accuracy': 0.9313858695652174,
    'negative_predictive_value': 0.9992266047950503}