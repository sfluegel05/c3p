"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids typically have a long aliphatic chain, with characteristic 
    functional groups such as 2-amino-1,3-diol backbones.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for sphingoid base - 2-amino-1,3-diol pattern
    # This is a minimal requirement and could vary greatly among different sphingoids
    sphingoid_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@@H](O)[CX4]")
    
    if mol.HasSubstructMatch(sphingoid_pattern):
        return True, "Contains sphingoid base pattern: 2-amino-1,3-diol"
    
    # Check for long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if mol.HasSubstructMatch(aliphatic_chain_pattern):
        return True, "Contains a long aliphatic chain, characteristic of sphingoids"
    
    return None, None  # Placeholder return if classification is ambiguous


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35785',
                          'name': 'sphingoid',
                          'definition': 'Sphinganine, its homologs and '
                                        'stereoisomers, and the hydroxy and '
                                        'unsaturated derivatives of these '
                                        'compounds.',
                          'parents': ['CHEBI:26739'],
                          'xrefs': ['LIPID_MAPS_class:LMSP01'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': None},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': None},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': None},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': None},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': None},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': None},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': None},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': None},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': None},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': None}],
    'sample_false_negatives': [   {   'smiles': 'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@@H](O)[C@H](C)N',
                                      'name': '1-deoxysphinga-6Z,9Z,12Z,15Z-tetraenine',
                                      'reason': None}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as "sphingoid" based on their '
                 "SMILES strings, we'll identify common structural features in "
                 'the provided examples. By examining the example SMILES, we '
                 'can formulate some core structural characteristics typical '
                 'of sphingoid compounds. These often include:\n'
                 '\n'
                 '1. A long aliphatic chain.\n'
                 '2. A backbone featuring hydroxyl groups and an amino group, '
                 'often in a 2-amino-1,3-diol configuration.\n'
                 '3. Often, one of the hydroxyl groups and the amino group are '
                 'adjacent, indicating a potential sphinganine or sphingosine '
                 'structure.\n'
                 '4. Potential additional functional groups, such as esters or '
                 'amides.\n'
                 '5. The presence of specific stereochemistry denoted as `@` '
                 'for chiral centers.\n'
                 '\n'
                 'Given these features, we can attempt to implement an '
                 'RDKit-based function to determine if a given SMILES string '
                 'matches the description of a "sphingoid" structure. The '
                 'search for characteristic substructures can be done using '
                 'SMARTS patterns. However, the complexity of sphingoid '
                 'structures and their diversity means that defining precise '
                 'patterns might be challenging, and thus if necessary, simply '
                 'return a generic identification message.\n'
                 '\n'
                 'Here is a Python code snippet utilizing RDKit to develop '
                 'such a classifier:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 94,
    'num_false_positives': 100,
    'num_true_negatives': 269,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.4845360824742268,
    'recall': 0.9894736842105263,
    'f1': 0.6505190311418685,
    'accuracy': 0.7823275862068966,
    'negative_predictive_value': 0.9962962962962963}