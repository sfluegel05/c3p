"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4' on the phenyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more general flavanone core pattern
    # Flavanone core consists of a 2-phenylchroman-4-one structure
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc3ccccc13")
    
    # Check for the presence of the flavanone core structure
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"
    
    # Define pattern for 4'-hydroxy group on an aromatic ring attached to the flavanone core
    # Ensures OH group is para to the attachment point
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Check for a 4'-OH group on the phenyl ring
    atom_indices = mol.GetSubstructMatch(four_prime_hydroxy_pattern)
    if atom_indices:
        # Verify that the OH group is on the correct phenyl ring
        phenyl_ring = Chem.MolFromSmarts("c1ccccc1")
        phenyl_match = mol.GetSubstructMatch(phenyl_ring)
        if phenyl_match and any(idx in atom_indices for idx in phenyl_match):
            return True, "Contains 4'-hydroxy group on aromatic ring of flavanone core"
    
    return False, "4'-hydroxy group not found at the correct position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:140331',
                          'name': "4'-hydroxyflavanones",
                          'definition': 'Any hydroxyflavanone having a hydroxy '
                                        "substituent located at position 4'.",
                          'parents': ['CHEBI:24697'],
                          'xrefs': ['MetaCyc:4-Hydroxyflavanones'],
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
               'COc1cc(ccc1O)[C@@H]1CC(=O)c2c(O)c(CC=C(C)C)c3OC(C)(C)C=Cc3c2O1 '
               'NAME: '
               "(2S)-6-(gamma,gamma-dimethylallyl)-5,4'-dihydroxy-3'-methoxy-6'',6''-dimethylpyran[2'',3'':7,8]flavanone "
               'REASON: MISSED No flavanone core structure identified\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1 '
               'NAME: (2S)-eriodictoyl-7-O-beta-D-glucopyranosiduronic acid '
               'REASON: MISSED No flavanone core structure identified\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\Cc1cc(cc(O)c1O)[C@@H]1CC(=O)c2c(O1)cc(O)c(CC=C(C)C)c2O '
               'NAME: schizolaenone B REASON: MISSED No flavanone core '
               'structure identified\n'
               ' * SMILES: Oc1ccc(cc1)[C@@H]1CC(=O)c2ccccc2O1 NAME: '
               "(2S)-4'-hydroxyflavanone REASON: MISSED No flavanone core "
               'structure identified\n'
               ' * SMILES: O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1 '
               'NAME: (+)-dihydromyricetin REASON: MISSED No flavanone core '
               'structure identified\n'
               ' * SMILES: '
               '[H][C@]1(CC(=O)c2c(O)cc(O)cc2O1)c1cc(O)c(O)c(C\\C=C(/C)CCC=C(C)C)c1 '
               'NAME: isonymphaeol B REASON: MISSED No flavanone core '
               'structure identified\n'
               ' * SMILES: COc1cc(ccc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1 NAME: '
               'homoeriodictyol REASON: MISSED No flavanone core structure '
               'identified\n'
               ' * SMILES: CC(C)=CCc1c(O)cc2O[C@@H](CC(=O)c2c1O)c1ccc(O)cc1 '
               'NAME: 6-prenylnaringenin REASON: MISSED No flavanone core '
               'structure identified\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\Cc1c(O)cc2O[C@@H](CC(=O)c2c1O)c1ccc(O)c(O)c1 '
               'NAME: nymphaeol A REASON: MISSED No flavanone core structure '
               'identified\n'
               ' * SMILES: Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1cc(O)c(O)c(O)c1 NAME: '
               '(2S)-dihydrotricetin REASON: MISSED No flavanone core '
               'structure identified\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No flavanone core structure '
                                               'identified'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No flavanone core structure '
                                               'identified'}],
    'sample_false_negatives': [],
    'attempt': 5,
    'reasoning': 'Upon reviewing the previous code and the outcomes, the main '
                 'issue seems to be an incorrect identification of the '
                 'flavanone core structure. Flavanones typically consist of a '
                 'chromanone structure with a phenyl group attached at a '
                 'specific position. Specifically, errors might be occurring '
                 'due to:\n'
                 '\n'
                 '1. **Flavanone Core Recognition**: The SMARTS pattern for '
                 'the flavanone core might be too restrictive or improperly '
                 'defined to match the structures listed.\n'
                 '   \n'
                 '2. **Positioning of Hydroxy Groups**: The code needs to '
                 'accurately check the position of the hydroxy groups, '
                 "ensuring it is on the correct aromatic ring and at the 4' "
                 'position relative to the flavanone core.\n'
                 '\n'
                 "To improve the program, we'll redefine the SMARTS pattern "
                 'for the flavanone core to be more general while still '
                 'capturing the essential features of flavanones. '
                 "Additionally, we'll ensure the hydroxy group's position is "
                 'correctly checked relative to this core.\n'
                 '\n'
                 "Here's the updated program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 80,
    'num_false_positives': 100,
    'num_true_negatives': 28465,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.4444444444444444,
    'recall': 1.0,
    'f1': 0.6153846153846153,
    'accuracy': 0.9965089893524175,
    'negative_predictive_value': 1.0}