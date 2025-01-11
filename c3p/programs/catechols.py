"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol or related moiety based on its SMILES string.
    Catechol is characterized by the presence of a 1,2-dihydroxybenzene moiety directly,
    or in a form that suggests such a phenolic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol or related moiety, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for catechol and related structures
    patterns = [
        Chem.MolFromSmarts("Oc1ccc(O)c1"),  # Core catechol (1,2-dihydroxybenzene)
        Chem.MolFromSmarts("Oc1c(O)cccc1"),  # To ensure overlooking by substituents
        Chem.MolFromSmarts("c1c(O)cccc1O")  # Other orientation for detecting phenolic OH
    ]

    # Check for presence of catechol or related patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains catechol or related moiety (detected)"

    return False, "Does not contain catechol or related moiety"

# Example usage
smiles_example = "Oc1ccc(O)c1"  # Catechol example
result, reason = is_catechols(smiles_example)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33566',
                          'name': 'catechols',
                          'definition': 'Any compound containing an o-diphenol '
                                        'component.',
                          'parents': ['CHEBI:33570'],
                          'xrefs': ['KEGG:C15571'],
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
               'O=C1OC[C@@H](C)[C@@H]([C@@H]1O)C2=C(O)C(O)=C(CO)C=C2 NAME: '
               'Calopin B REASON: MISSED Does not contain catechol moiety\n'
               ' * SMILES: COc1cc(CCc2ccc(O)c(O)c2)cc(O)c1O NAME: dendrocandin '
               'E REASON: MISSED Does not contain catechol moiety\n'
               ' * SMILES: Oc1cccc(-c2ccc(Cl)cc2)c1O NAME: '
               "4'-chlorobiphenyl-2,3-diol REASON: MISSED Does not contain "
               'catechol moiety\n'
               ' * SMILES: C1(C2=CC=CC=C2O)=CC=CC(=C1O)O NAME: '
               "biphenyl-2,2',3-triol REASON: MISSED Does not contain catechol "
               'moiety\n'
               ' * SMILES: OC1=CC(=O)OC(\\C=C\\C2=CC=C(O)C(O)=C2)=C1 NAME: '
               'hispidin REASON: MISSED Does not contain catechol moiety\n'
               ' * SMILES: '
               'C1(C2=C(C(=CC(=C2)C(O)=O)O)O)=C(C(=CC(=C1)C(=O)O)OC)O NAME: '
               "5,6,6'-trihydroxy-5'-methoxy[biphenyl]-3,3'-dicarboxylic acid "
               'REASON: MISSED Does not contain catechol moiety\n'
               ' * SMILES: '
               'C12=CC(C3=C(C=C(C(=C3C)O)O)[C@@]1(CC[C@@]4([C@@]2(CC[C@@]5([C@]4(C[C@@](CC5)(C(O)=O)C)[H])C)C)C)C)=O '
               'NAME: wilforol A REASON: MISSED Does not contain catechol '
               'moiety\n'
               ' * SMILES: '
               'O=C1OC(=CC(=C1C(C=2C(=O)OC(/C=C/C3=CC(O)=C(O)C=C3)=CC2O)CCC(=O)OC)O)/C=C/C4=CC(O)=C(O)C=C4 '
               'NAME: Phaeolschidin D REASON: MISSED Does not contain catechol '
               'moiety\n'
               ' * SMILES: Oc1ccc(F)cc1O NAME: 4-Fluorocatechol REASON: MISSED '
               'Does not contain catechol moiety\n'
               ' * SMILES: [H]C(=O)Cc1ccc(O)c(O)c1 NAME: '
               '3,4-dihydroxyphenylacetaldehyde REASON: MISSED Does not '
               'contain catechol moiety\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Does not contain catechol or '
                                               'related moiety'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@H]1C=CC=C[C@H]1O',
                                      'name': 'cis-cyclohexa-3,5-diene-1,2-diol',
                                      'reason': 'Does not contain catechol or '
                                                'related moiety'},
                                  {   'smiles': 'O[C@H]1C=CC=C[C@@H]1O',
                                      'name': '(1S,2S)-cyclohexa-3,5-diene-1,2-diol',
                                      'reason': 'Does not contain catechol or '
                                                'related moiety'},
                                  {   'smiles': 'Oc1ccc(O)c(c1)S(O)(=O)=O',
                                      'name': '2,5-dihydroxybenzenesulfonic '
                                              'acid',
                                      'reason': 'Does not contain catechol or '
                                                'related moiety'}],
    'attempt': 1,
    'reasoning': 'The previous implementation provided a reasonable start by '
                 'correctly identifying the presence of a catechol moiety, '
                 'defined as 1,2-dihydroxybenzene. However, based on the '
                 'outcomes described, the false negatives suggest that the '
                 "pattern recognition for catechols could be improved. Let's "
                 'consider the following improvements to address these '
                 'issues:\n'
                 '\n'
                 '1. **Pattern Specificity vs. Generality**: The SMARTS '
                 'pattern `Oc1ccc(O)c1` may need refining or additional '
                 'patterns to adequately capture variations in the catechol '
                 'moiety, such as meta-substitution or fused rings carrying '
                 'the catechol moiety.\n'
                 '\n'
                 '2. **Structural Variability**: Molecules may contain '
                 'catechol-related structures like substituted catechols or '
                 'fused ring systems that carry the catechol pattern. Ensure '
                 'the SMARTS pattern captures these possibilities.\n'
                 '\n'
                 '3. **Complex Structures Support**: Compounds like '
                 '`dendrocandin E` involve complex structures. Additional or '
                 'extended checks around the primary catechol determinations '
                 'may be necessary due to extended conjugation, linkages, '
                 'etc.\n'
                 '\n'
                 '4. **Update SMARTS Patterns**: SMARTS might need to be '
                 'adjusted by incorporating alternate representations, '
                 'ensuring inclusivity of inherent/subtle structures which '
                 'maintain catechol-like activity. This might mean adding more '
                 'conditions or structures that specifically avoid excluded '
                 'groups.\n'
                 '\n'
                 'Based on this reasoning, I’m going to update the program to '
                 'improve its accuracy by expanding the SMARTS patterns used '
                 'for matching catechols and offering comments and checks for '
                 'each step:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 322,
    'num_false_positives': 100,
    'num_true_negatives': 909,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.7630331753554502,
    'recall': 0.9907692307692307,
    'f1': 0.8621151271753682,
    'accuracy': 0.9227886056971514,
    'negative_predictive_value': 0.9967105263157895}