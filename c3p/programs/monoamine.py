"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine typically has an amino group linked to an aromatic ring by a two-carbon bridge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aromatic ring pattern
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")  # benzene as a simple example
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"
    
    # Look for two-carbon alkyl chain connected to an amino group
    amine_pattern = Chem.MolFromSmarts("NCC")  # simple primary or secondary amine pattern
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No two-carbon chain linked to an amino group found"

    return True, "Contains an aromatic ring and an amino group linked via a two-carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63534',
                          'name': 'monoamine',
                          'definition': 'An aralylamino compound which '
                                        'contains one amino group connected to '
                                        'an aromatic ring by a two-carbon '
                                        'chain. Monoamines are derived from '
                                        'aromatic amino acids like '
                                        'phenylalanine, tyrosine, tryptophan, '
                                        'and the thyroid hormones by the '
                                        'action of aromatic amino acid '
                                        'decarboxylase enzymes.',
                          'parents': ['CHEBI:64365'],
                          'xrefs': [   'PMID:21822758',
                                       'PMID:21993877',
                                       'PMID:22005599',
                                       'PMID:22082101',
                                       'PMID:22153577',
                                       'PMID:22213370',
                                       'PMID:22218931',
                                       'PMID:22342987',
                                       'PMID:22371656'],
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
                                     'reason': 'No two-carbon chain linked to '
                                               'an amino group found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No two-carbon chain linked to '
                                               'an amino group found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'No aromatic ring found'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'No two-carbon chain linked to '
                                               'an amino group found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as monoamines based on their '
                 'SMILES string, we need to identify key structural features '
                 'of monoamines. Monoamines are organic compounds '
                 'characterized by having one amino group connected to an '
                 'aromatic ring via a two-carbon alkyl chain (ethylene '
                 'bridge).\n'
                 '\n'
                 'Key features to look for in the structure include:\n'
                 '1. An aromatic ring system (e.g., benzene).\n'
                 '2. An amino group (-NH2 or a substituted version like -NH- '
                 'or -N<).\n'
                 '3. A two-carbon chain linking the amino group to the '
                 'aromatic ring.\n'
                 '\n'
                 "We will use RDKit's SMILES and SMARTS parsing capabilities "
                 'to identify these patterns.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 33,
    'num_false_positives': 100,
    'num_true_negatives': 232,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.24812030075187969,
    'recall': 1.0,
    'f1': 0.3975903614457831,
    'accuracy': 0.726027397260274,
    'negative_predictive_value': 1.0}