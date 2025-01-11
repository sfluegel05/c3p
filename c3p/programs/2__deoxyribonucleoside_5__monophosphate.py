"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    These molecules have a 2'-deoxyribose sugar, a phosphate group at the 5' position,
    and a nucleobase attached to the 1' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")
    
    # Check for 2'-deoxyribose moiety (allowing variable stereochemistry)
    deoxyribose_pattern = Chem.MolFromSmarts('[C@H]1C[C@@H](O)[C@H](CO)O1')
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return (False, "No 2'-deoxyribose moiety found")

    # Check for a phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts('COP(=O)(O)O')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return (False, "No 5'-phosphate group found")
        
    # Check for nucleobase attached to the sugar via a variable pattern
    nucleobase_pattern = Chem.MolFromSmarts('n')
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return (False, "No nucleobase found (purine/pyrimidine) attached to the sugar")
    
    return (True, "Contains 2'-deoxyribose, 5'-phosphate group, and a nucleobase")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18241',
                          'name': "2'-deoxyribonucleoside 5'-monophosphate",
                          'definition': "A 2'-deoxyribonucleoside "
                                        'monophosphate compound with the '
                                        "phosphate group in the 5'-position.",
                          'parents': ['CHEBI:19257', 'CHEBI:37016'],
                          'xrefs': ['KEGG:C00676', 'KEGG:C03607'],
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
               'C12C(C3C1N(C(N=C3N)=O)[C@@H]4O[C@@H]([C@H](C4)O)COP(=O)(O)O)C(=NC(N2[C@@H]5O[C@@H]([C@H](C5)O)COP(=O)(O)O)=O)N '
               "NAME: 2'-deoxycytidine 5'-monophosphate dimer REASON: MISSED "
               "No 2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'Cc1cn([C@@H]2C[C@@H](O)[C@H](COP(O)(O)=O)O2)c(=O)[nH]c1=O '
               'NAME: 1-(2-deoxy-5-O-phosphono-beta-L-ribofuranosyl)thymine '
               "REASON: MISSED No 2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(O)=O)O1 NAME: '
               "2'-deoxyadenosine 5'-monophosphate REASON: MISSED No "
               "2'-deoxyribose moiety found\n"
               ' * SMILES: Nc1ccn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)n1 '
               "NAME: 2'-deoxycytosine 5'-monophosphate REASON: MISSED No "
               "2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'Nc1nc2n([C@H]3C[C@H](O)[C@@H](COP(O)(O)=O)O3)c(=O)[nH]c2c(=O)[nH]1 '
               "NAME: 8-oxo-dGMP REASON: MISSED No 2'-deoxyribose moiety "
               'found\n'
               ' * SMILES: '
               'O[C@H]1C[C@@H](O[C@@H]1COP(O)(O)=O)n1ccc(=O)[nH]c1=O NAME: '
               "dUMP REASON: MISSED No 2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'NC1=NC=NC2=C1N=NN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](O)C3 NAME: '
               "8-aza-2'-deoxyadenosine 5'-monophosphate REASON: MISSED No "
               "2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'C1(=C(O)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O)O '
               "NAME: 5,6-dihydroxy-2'-deoxyuridine 5'-monophosphate REASON: "
               "MISSED No 2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'C1C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O NAME: '
               "5,6-dihydrothymidine 5'-monophosphate REASON: MISSED No "
               "2'-deoxyribose moiety found\n"
               ' * SMILES: '
               'C1=C(CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)C(NC(N1[C@@H]3O[C@H](COP(=O)(O)O)[C@H](C3)O)=O)=O '
               'NAME: (beta-D-glucopyranosyloxymethyl)deoxyuridine '
               "5'-monophosphate REASON: MISSED No 2'-deoxyribose moiety "
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': "No 2'-deoxyribose moiety "
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C12C(C3C1N(C(N=C3N)=O)[C@@H]4O[C@@H]([C@H](C4)O)COP(=O)(O)O)C(=NC(N2[C@@H]5O[C@@H]([C@H](C5)O)COP(=O)(O)O)=O)N',
                                      'name': "2'-deoxycytidine "
                                              "5'-monophosphate dimer",
                                      'reason': 'No nucleobase found '
                                                '(purine/pyrimidine) attached '
                                                'to the sugar'},
                                  {   'smiles': 'C1C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O',
                                      'name': '5,6-dihydrothymidine '
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase found '
                                                '(purine/pyrimidine) attached '
                                                'to the sugar'},
                                  {   'smiles': 'C1(C(C(=NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)N)O)O',
                                      'name': "2'-deoxy-5,6-dihydroxy-5,6-dihydrocytidine "
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase found '
                                                '(purine/pyrimidine) attached '
                                                'to the sugar'},
                                  {   'smiles': 'C12C(C3(C1N(C(NC3=O)=O)[C@@H]4O[C@@H]([C@H](C4)O)COP(=O)(O)O)C)C(=NC(N2[C@@H]5O[C@@H]([C@H](C5)O)COP(=O)(O)O)=O)N',
                                      'name': 'dCMP-dTMP',
                                      'reason': 'No nucleobase found '
                                                '(purine/pyrimidine) attached '
                                                'to the sugar'},
                                  {   'smiles': 'C12N(C(NC(C1(C)C3(C2N(C(NC3=O)=O)[C@@H]4O[C@@H]([C@H](C4)O)COP(=O)(O)O)C)=O)=O)[C@@H]5O[C@@H]([C@H](C5)O)COP(=O)(O)O',
                                      'name': "thymidine 5'-monophosphate "
                                              'dimer',
                                      'reason': 'No nucleobase found '
                                                '(purine/pyrimidine) attached '
                                                'to the sugar'}],
    'attempt': 1,
    'reasoning': "To address the failure in the previous attempt, let's "
                 'analyze and adjust the conditions for identifying a '
                 "2'-deoxyribonucleoside 5'-monophosphate:\n"
                 '\n'
                 '### Key structural features:\n'
                 "1. **2'-deoxyribose sugar**: Unlike ribose, this sugar lacks "
                 "an -OH group at the C2' position. It's important to ensure "
                 "we're flexible with stereochemistry, given that both beta "
                 'and alpha configurations might exist, and different '
                 'orientations in the SMILES.\n'
                 "2. **Phosphate group at 5' position**: Appears as "
                 '`COP(=O)(O)O` in SMILES.\n'
                 "3. **Nucleobase attached to the sugar base at C1'**: The "
                 'base can be different, including purine or pyrimidine '
                 'manufacturers. It might be directly or indirectly attached '
                 'due to various substitutions or ring alterations.\n'
                 '\n'
                 '### Corrective actions:\n'
                 '- **Broaden stereochemistry flexibility**: The previous '
                 'pattern was too restrictive concerning stereochemistry.\n'
                 '- **Adjust nucleobase detection**: We need a more flexible '
                 'representation for nucleobases since they might not exactly '
                 'match the given pattern due to modifications.\n'
                 '\n'
                 '### Improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 20,
    'num_false_positives': 100,
    'num_true_negatives': 4902,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.16666666666666666,
    'recall': 0.8,
    'f1': 0.27586206896551724,
    'accuracy': 0.9791127909289835,
    'negative_predictive_value': 0.9989810474831873}