"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA derivative based on its SMILES string.

    A 2-enoyl-CoA derivative features a CoA moiety (Coenzyme A), typically recognized by
    an adenosine-linked phosphate chain, linked via thioester to a fatty acyl chain with an enoyl group (C=C).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 2-enoyl-CoA derivative with a reason, False otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified CoA part pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O")
    if not coa_pattern or not mol.HasSubstructMatch(coa_pattern):
        return False, "No or incomplete CoA moiety found"
    
    # Enoyl group pattern (C=C double bond in acyl chain)
    enoyl_pattern = Chem.MolFromSmarts("C=C")
    if not enoyl_pattern or not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No enoyl group found"

    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not thioester_pattern or not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Pass all checks
    return True, "Contains CoA moiety, enoyl group, and thioester linkage"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:19573',
                          'name': '2-enoyl-CoA',
                          'definition': 'An unsaturated fatty acyl-CoA in '
                                        'which the S-acyl group contains a '
                                        'double bond between positions 2 and '
                                        '3.',
                          'parents': ['CHEBI:51006'],
                          'xrefs': ['KEGG:C00605', 'MetaCyc:D2-ENOYL-COA'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.HasSubstructMatch(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params=True)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool recursionPossible=True, bool useChirality=False, '
               'bool useQueryQueryMatches=False)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'bool recursionPossible=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False)\n'
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
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No or incomplete CoA moiety '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'In the previous classification attempt for 2-enoyl-CoA, the '
                 'function aimed to identify molecules with a CoA moiety, an '
                 'enoyl group, and a thioester linkage. However, the function '
                 'might have failed due to several reasons:\n'
                 '\n'
                 '1. **CoA Pattern Complexity:** The defined SMARTS pattern '
                 'for identifying CoA moiety may have been too specific or not '
                 'accurately representing all variations of CoA structures '
                 'seen in 2-enoyl-CoA derivatives. This could result in the '
                 'function failing to identify valid CoA structures.\n'
                 '\n'
                 '2. **Matching Errors:** Errors in substructure matching, '
                 'such as the one indicated by the log, suggest that SMARTS '
                 'patterns might not be valid or constructed correctly. There '
                 'is a need for careful validation and possible simplification '
                 'of SMARTS patterns to correctly match substructures in '
                 'molecules.\n'
                 '\n'
                 '3. **Lack of Flexibility for Structural Variability:** The '
                 'function might not account for possible stereochemical and '
                 'configurational variations present in the structures, '
                 'causing valid 2-enoyl-CoA derivatives to be missed or '
                 'incorrectly rejected.\n'
                 '\n'
                 'To improve the classification, consider simplifying and '
                 'validating the SMARTS patterns to ensure they match typical '
                 'components without being too restrictive, while also '
                 'accommodating structural variation:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 43,
    'num_false_positives': 100,
    'num_true_negatives': 21842,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3006993006993007,
    'recall': 1.0,
    'f1': 0.46236559139784944,
    'accuracy': 0.9954514441664771,
    'negative_predictive_value': 1.0}