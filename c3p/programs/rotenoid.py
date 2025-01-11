"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are characterized by a polycyclic structure often with chromene or isoflavone-like backbones with methoxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for rotenoid substructures
    # A proposed pattern might include a chromene unit and a methoxy substituent
    chromene_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")
    methoxy_pattern = Chem.MolFromSmarts("CO")

    # Check for the presence of a chromene-like backbone
    if not mol.HasSubstructMatch(chromene_pattern):
        return False, "No chromene-like backbone found"

    # Check for methoxy substituents
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 1:
        return False, "No methoxy groups found"
    
    # Further checks can be done to match additional known rotenoid features (optional)
    # e.g., additional oxygen heterocycles, glycoside linkages, etc.
    
    return True, "Contains chromene-like backbone with methoxy group(s)"

# Example use
# an example SMILES for a rotenoid could be validated here:
# result, reason = is_rotenoid("SMILES HERE")
# print(f"Is rotenoid: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71543',
                          'name': 'rotenoid',
                          'definition': 'Members of the class of '
                                        'tetrahydrochromenochromene that '
                                        'consists of a cis-fused '
                                        'tetrahydrochromeno[3,4-b]chromene '
                                        'skeleton and its substituted '
                                        'derivatives. The term was originally '
                                        'restricted to natural products, but '
                                        'is now also used to describe '
                                        'semi-synthetic and fully synthetic '
                                        'compounds.',
                          'parents': ['CHEBI:72544', 'CHEBI:72579'],
                          'xrefs': ['Wikipedia:Rotenoids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No chromene-like backbone '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No chromene-like backbone '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12COC3=C(C=C(OC)C(OC)=C3)[C@]1([H])C(=O)C1=C(O2)C2=C(OC(C)(C)C=C2)C=C1',
                                      'name': 'deguelin',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O1C2C(O)(C(O)C3=C1C=4CC(OC4C=C3)C(CO[C@@H]5OC([C@@H](O)[C@H](O)C5O)CO)=C)C=6C(OC2)=CC(OC)=C(OC)C6',
                                      'name': '12-Dihydrodalbinol O-glucoside',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O1C2C(C(=O)C3=C1C=4CC(OC4C=C3)C(O)(CO)C)C=5C(OC2)=CC(OC)=C(OC)C5',
                                      'name': 'Amorphigenol',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': '[H][C@@]12COc3cc(OC)c(OC)cc3[C@]1(O)C(=O)c1ccc3OC(C)(C)C=Cc3c1O2',
                                      'name': 'tephrosin',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O1C2=C(C=3C(OC2OC)=CC=CC3)C(=O)C4=C1C=C(O)C(=C4O)C',
                                      'name': 'Boeravinone A',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': '[H][C@@]12COc3cc(OC)c(OC)cc3[C@]1(O)C(=O)c1c(O)cc3OC(C)(C)C(O)C(OC)c3c1O2',
                                      'name': "4',5'-dihydro-11,5'-dihydroxy-4'-methoxytephrosin",
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O(C(C1OC2=C(C1)C=3OC4C(C(=O)C3C=C2)C=5C(OC4)=CC(OC)=C(OC)C5)(CO)C)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO[C@@H]7OC[C@H](O)C(O)C7O',
                                      'name': 'Amorphigenol O-vicianoside',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O(C(C1OC2=C(C1)C=3OC4C(C(=O)C3C=C2)C=5C(OC4)=CC(OC)=C(OC)C5)(CO)C)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO',
                                      'name': 'Amorphigenol O-glucoside',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O1C2C(O)(C=3C(OC2)=CC(OC)=C(OC)C3)C(=O)C4=C1C=5CC(OC5C=C4)C(C)=C',
                                      'name': '13-Hydroxy-16,17-dimethoxy-6-prop-1-en-2-yl-2,7,20-trioxapentacyclo[11.8.0.03,11.04,8.014,19]henicosa-3(11),4(8),9,14,16,18-hexaen-12-one',
                                      'reason': 'No chromene-like backbone '
                                                'found'},
                                  {   'smiles': 'O1C2C(C(=O)C3=C1C=4CC(OC4C=C3)C(CO)C)C=5C(OC2)=CC(OC)=C(OC)C5',
                                      'name': 'Dihydroamorphigenin',
                                      'reason': 'No chromene-like backbone '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify whether a chemical entity is a rotenoid based on '
                 'its SMILES string, we should first understand the structural '
                 'characteristics that define such compounds. Rotenoids are '
                 'naturally occurring plant compounds, known for having a '
                 'complex polycyclic structure with oxygen heterocycles. '
                 'Typically, they have a core structure derived from a '
                 'polycyclic system that may include a chromene or '
                 'isoflavone-like backbone, often functionalized with one or '
                 'more methoxy groups, and sometimes additional rings or '
                 'glycosidic linkages.\n'
                 '\n'
                 'Our task is to create a function `is_rotenoid` that uses '
                 'SMARTS (SMILES Arbitrary Target Specification) patterns to '
                 'identify these structural features in a given SMILES string. '
                 'The function will leverage RDKit to match these patterns.\n'
                 '\n'
                 "Here's what our code will do:\n"
                 '\n'
                 '1. Use RDKit to parse the SMILES string to a molecule '
                 'representation.\n'
                 '2. Define SMARTS patterns that capture key structural '
                 'features of rotenoids - for instance, a chromene backbone '
                 'and methoxy groups.\n'
                 '3. Use RDKit to check if a molecule contains these '
                 'substructures. \n'
                 '4. If matches for these essential features are found, '
                 'classify the molecule as a rotenoid. Otherwise, classify it '
                 'as non-rotenoid.\n'
                 '\n'
                 "Let's proceed with the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 8,
    'num_false_positives': 100,
    'num_true_negatives': 30014,
    'num_false_negatives': 28,
    'num_negatives': None,
    'precision': 0.07407407407407407,
    'recall': 0.2222222222222222,
    'f1': 0.1111111111111111,
    'accuracy': 0.9957545605306799,
    'negative_predictive_value': 0.9990679715065575}