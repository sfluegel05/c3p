"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"
    
    # Look for a primary alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[CX4;H2][OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No primary alcohol (C-OH) group found"

    return True, "Contains an aromatic ring and a primary alcohol group"

# Testing example
# is_aromatic_primary_alcohol("C1=CC(=CN=C1)CO")  # 3-pyridinemethanol, expected True.


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33857',
                          'name': 'aromatic primary alcohol',
                          'definition': 'Any primary alcohol in which the '
                                        'alcoholic hydroxy group is attached '
                                        'to a carbon which is itself bonded to '
                                        'an aromatic ring.',
                          'parents': ['CHEBI:15734', 'CHEBI:33854'],
                          'xrefs': ['KEGG:C03485'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No primary alcohol (C-OH) '
                                               'group found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No primary alcohol (C-OH) '
                                               'group found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No primary alcohol (C-OH) '
                                               'group found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No primary alcohol (C-OH) '
                                               'group found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'No aromatic ring found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as aromatic primary alcohols '
                 'using their SMILES string, we need a strategy to identify '
                 'key structural features. An aromatic primary alcohol '
                 'typically includes:\n'
                 '\n'
                 '1. An aromatic ring, which could be a benzene ring or a '
                 'heteroaromatic ring containing atoms like nitrogen, oxygen, '
                 'or sulfur.\n'
                 '2. A primary alcohol group, which consists of a hydroxyl '
                 'group (OH) directly attached to the primary carbon that is '
                 'not part of the aromatic ring.\n'
                 '\n'
                 'For this task, we need to use RDKit to help us identify '
                 "these features through SMARTS patterns. Here's how the "
                 'process should be laid out:\n'
                 '\n'
                 "1. **Identify the aromatic ring**: We'll use a SMARTS "
                 'pattern like `a` to check for any aromatic atoms since '
                 "RDKit's SMILES parser can classify atoms as aromatic.\n"
                 "2. **Identify the primary alcohol group**: We'll search for "
                 'a SMARTS pattern like `CO` where `C` is a primary carbon sp3 '
                 'hybridized, and the `O` in `OH` is connected.\n'
                 '\n'
                 'The algorithm will check for these patterns in the given '
                 'SMILES and return whether the molecule fits the criteria of '
                 'an aromatic primary alcohol.\n'
                 '\n'
                 "Here's the step-by-step code implementation with additional "
                 'comments:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 36,
    'num_false_positives': 100,
    'num_true_negatives': 544,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2647058823529412,
    'recall': 1.0,
    'f1': 0.4186046511627907,
    'accuracy': 0.8529411764705882,
    'negative_predictive_value': 1.0}