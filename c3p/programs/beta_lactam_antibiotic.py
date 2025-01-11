"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic typically contains a beta-lactam ring: a 4-membered ring with one nitrogen atom.
    We also consider scaffolds characteristic of known antibiotics, such as penicillins, cephalosporins, and carbapenems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for beta-lactam antibiotics
    beta_lactam_patterns = [
        Chem.MolFromSmarts("C1CNC(=O)1"),  # General beta-lactam ring
        Chem.MolFromSmarts("C1[C@@H](N1)C(=O)[C@H]"),  # Stereochemically specific beta-lactam
        
        # Penicillin core: thiazolidine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC[C@](N1)C(=O)"),  # Thiazolidine ring with chirality
        
        # Cephalosporin core: dihydrothiazine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC2CCN1C2=O"),  # Dihydrothiazine ring with chirality
        
        # Carbapenem core: unsaturated beta-lactam with sulfur
        Chem.MolFromSmarts("C1C(=C)N(C1=O)C(=O)"),
    ]
    
    # Check if any of the beta-lactam patterns match
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a beta-lactam ring in a known antibiotic scaffold"
    
    return False, "No beta-lactam ring found or not in a typical antibiotic context"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27933',
                          'name': 'beta-lactam antibiotic',
                          'definition': 'An organonitrogen heterocyclic '
                                        'antibiotic that contains a '
                                        'beta-lactam ring.',
                          'parents': ['CHEBI:25558', 'CHEBI:35627'],
                          'xrefs': [   'KEGG:C03438',
                                       'PMID:19254642',
                                       'PMID:22594007',
                                       'Wikipedia:Beta-lactam_antibiotic'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: invalid syntax. Perhaps you forgot a comma? (<string>, '
               'line 32)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No beta-lactam ring found or '
                                               'not in a typical antibiotic '
                                               'context'}],
    'sample_false_negatives': [   {   'smiles': '[H]C(N[C@@H](C([O-])=O)C(C)(C)S)=C1N=C(Cc2ccccc2)OC1=O',
                                      'name': 'benzylpenicillenate',
                                      'reason': 'No beta-lactam ring found or '
                                                'not in a typical antibiotic '
                                                'context'}],
    'attempt': 3,
    'reasoning': 'To improve the extraction of beta-lactam antibiotics based '
                 'on their SMILES strings, it is important to address the '
                 'mistakes in the previous code and enhance it by considering '
                 'specific patterns found in known antibiotics composed of a '
                 'beta-lactam ring. The errors possibly arose from incorrect '
                 'use of SMARTS patterns and syntax issues.\n'
                 '\n'
                 '### Key Steps to Consider:\n'
                 '1. **Correcting Syntax Error**: Ensure proper syntax, '
                 'especially missing commas between elements in an array like '
                 '`beta_lactam_patterns`.\n'
                 '2. **SMARTS Pattern Definition**:\n'
                 '   - **General Beta-Lactam Core**: 4-membered ring with a '
                 'nitrogen and a carbonyl group is essential.\n'
                 '   - **Known Antibiotic Scaffolds**: Patterns specific to '
                 'penicillins, cephalosporins, and carbapenems can aid in '
                 'further identification.\n'
                 '3. **Pattern Matching on Molecule**: Use RDKit to handle the '
                 'molecule parsing and SMARTS pattern matching correctly.\n'
                 '4. **Explanatory Messaging**: Provide clear feedback on '
                 'classification outcome.\n'
                 '\n'
                 '### Revised Program:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 253,
    'num_false_positives': 100,
    'num_true_negatives': 72765,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.71671388101983,
    'recall': 0.9960629921259843,
    'f1': 0.8336079077429984,
    'accuracy': 0.9986186900805536,
    'negative_predictive_value': 0.9999862573179782}