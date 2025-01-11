"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a structure with a sugar-like hydroxy pattern
    and one or more phosphate groups esterified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefining a pattern for sugar-like hydroxy structures
    # This accounts for both furanose and pyranose ring sugars
    sugar_ring_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1O")  # Typical hexose as an example
    # Linear sugar parts (e.g., glyceraldehyde phosphate)
    sugar_linear_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # Refine patterns for phosphate groups (covering monophosphate, diphosphate, etc.)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O])[O]")  # Phosphate group

    # Check for presence of sugar structure
    has_sugar_structure = mol.HasSubstructMatch(sugar_ring_pattern) or mol.HasSubstructMatch(sugar_linear_pattern)
    if not has_sugar_structure:
        return False, "No sugar-like structure found"

    # Check for presence of phosphate group
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group esterified found"

    # If both structures are present, it's a phospho sugar
    return True, "Contains a sugar-like structure with phosphate esterification"

# Example metadata for the classification
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:??',
        'name': 'phospho_sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    },
    'message': None,
    'success': True
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33447',
                          'name': 'phospho sugar',
                          'definition': 'Any monosaccharide containing an '
                                        'alcoholic hydroxy group esterified '
                                        'with phosphoric acid.',
                          'parents': ['CHEBI:26816', 'CHEBI:63367'],
                          'xrefs': ['KEGG:C00934', 'PMID:18186488'],
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
               'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C[C@H]1O '
               'NAME: cordycepin triphosphate REASON: MISSED No polyhydroxy '
               'structure typical of monosaccharides found\n'
               ' * SMILES: '
               'N[C@@H](CO)C(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12 '
               'NAME: L-seryl-AMP REASON: MISSED No polyhydroxy structure '
               'typical of monosaccharides found\n'
               ' * SMILES: '
               'O[C@H]1C[C@@H](O[C@@H]1COP(O)(=O)OP(O)(=O)OP(O)(O)=O)N1C=CC(=O)NC1=O '
               'NAME: dUTP REASON: MISSED No polyhydroxy structure typical of '
               'monosaccharides found\n'
               ' * SMILES: '
               'CCCC(=O)Nc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O '
               'NAME: N(6)-butyryl-cAMP REASON: MISSED No polyhydroxy '
               'structure typical of monosaccharides found\n'
               ' * SMILES: '
               'Cc1cc2[N]c3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)c2cc1C '
               'NAME: FMNH(.) REASON: MISSED No polyhydroxy structure typical '
               'of monosaccharides found\n'
               ' * SMILES: '
               'Cc1nc2c(N)ncnc2n1[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O '
               "NAME: 8-methyladenosine 5'-monophosphate REASON: MISSED No "
               'polyhydroxy structure typical of monosaccharides found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1cc(CNCC(O)=O)c(=O)[nH]c1=O '
               "NAME: 5-carboxymethylaminomethyluridine 5'-monophosphate "
               'REASON: MISSED No polyhydroxy structure typical of '
               'monosaccharides found\n'
               ' * SMILES: OC1C[C@H](O)[C@@H](COP(O)(O)=O)O1 NAME: '
               '2-deoxy-D-ribofuranose 5-phosphate REASON: MISSED No '
               'polyhydroxy structure typical of monosaccharides found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O '
               "NAME: uridine 5'-monophosphate REASON: MISSED No polyhydroxy "
               'structure typical of monosaccharides found\n'
               ' * SMILES: '
               'C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O '
               "NAME: 8-formyl-8-demethylriboflavin 5'-phosphate REASON: "
               'MISSED No polyhydroxy structure typical of monosaccharides '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No phosphate group esterified '
                                               'found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No sugar-like structure found'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)C(O)=O)[C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](OC(=O)C)[C@H](O)CO)C(O)=O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-6-[(1S,2R)-2-[(2S,4S,5R,6R)-5-acetamido-6-[(1R,2R)-1-acetyloxy-2,3-dihydroxypropyl]-2-carboxy-4-hydroxyoxan-2-yl]oxy-1,3-dihydroxypropyl]-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No phosphate group esterified '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C[C@H]1O',
                                      'name': 'cordycepin triphosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': '[H]C(=O)C(O)COP(O)(O)=O',
                                      'name': 'glyceraldehyde 3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCOCC(=O)COP(O)(O)=O',
                                      'name': '1-heptadecylglycerone '
                                              '3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)OCC(=O)COP(O)(O)=O',
                                      'name': '1-arachidonoylglycerone '
                                              '3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(O)=O)C(=O)[C@H]1O',
                                      'name': "3'-dehydro-AMP",
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'OCC(O)C(=O)COP(O)(O)=O',
                                      'name': 'erythrulose 1-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOCC(=O)COP(O)(O)=O',
                                      'name': '1-palmitylglycerone 3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': '[H]C(=O)[C@H](O)COP(O)(O)=O',
                                      'name': 'D-glyceraldehyde 3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'CC\\C=C/C\\C=C/C\\C=C/CCCCCCCCOCC(=O)COP(O)(O)=O',
                                      'name': '1-(9Z,12Z,15Z)-octadecatrienylglycerone '
                                              '3-phosphate',
                                      'reason': 'No sugar-like structure '
                                                'found'},
                                  {   'smiles': 'P(OCC(=O)COC(=O)CCCCCCCCC)(O)(O)=O',
                                      'name': 'DHAP(10:0)',
                                      'reason': 'No sugar-like structure '
                                                'found'}],
    'attempt': 2,
    'reasoning': 'To determine why the previous attempt at classifying phospho '
                 "sugars failed, let's summarize a few important points:\n"
                 '\n'
                 '1. **Previous Pattern Definitions**: The earlier code '
                 'defined patterns to identify polyhydroxy structures and '
                 'phosphate groups. The polyhydroxy pattern focused on '
                 'structures resembling sugars, and the phosphate pattern '
                 'looked for phosphate ester groups.\n'
                 '\n'
                 '2. **Failure Analysis**: The outcomes indicate all cases '
                 'were misclassified because the so-called "polyhydroxy '
                 'backbone" is too broad and did not capture the diversity of '
                 'possible monosaccharide structures in the given examples.\n'
                 '\n'
                 '3. **Improvement Plan**: \n'
                 '   - **Redefine Polyhydroxy Pattern**: The SMARTS pattern '
                 'should closely align with known sugar ring systems and '
                 'linear chains typical of phospho sugars.\n'
                 '   - **Extend Phosphate Pattern**: Ensure the phosphate '
                 'group pattern accounts for the potential complexity and '
                 'variation in phospho sugared forms.\n'
                 '   - **Energy-efficient Search**: Use more specific patterns '
                 'to minimize false positives.\n'
                 '\n'
                 'The following code block will address these issues by '
                 'refining the pattern definitions and checks:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 434,
    'num_false_positives': 100,
    'num_true_negatives': 1879,
    'num_false_negatives': 26,
    'num_negatives': None,
    'precision': 0.8127340823970037,
    'recall': 0.9434782608695652,
    'f1': 0.8732394366197181,
    'accuracy': 0.948339483394834,
    'negative_predictive_value': 0.9863517060367454}