"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 15 carbons (sesquiterpene)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 15:
        return False, f"Number of carbons is {c_count}, not 15 (not a sesquiterpene)"

    # Check for 10-membered ring (germacrane skeleton)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if 10 not in ring_sizes:
        return False, "No 10-membered ring found (germacrane skeleton missing)"

    # Check for lactone group (cyclic ester within a ring)
    lactone_smarts = "[C;R](=O)[O;R][C;R]"  # Ester group where atoms are in a ring
    lactone_mol = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_mol):
        return False, "No lactone ring (cyclic ester) found"

    return True, "Molecule is a germacranolide (15 carbons, 10-membered ring, and lactone group)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73011',
                          'name': 'germacranolide',
                          'definition': 'A sesquiterpene lactone based on '
                                        'germacrane skeleton.',
                          'parents': ['CHEBI:37667'],
                          'xrefs': ['Wikipedia:Germacranolide'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Number of carbons is 6, not 15 '
                                               '(not a sesquiterpene)'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Number of carbons is 13, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Number of carbons is 30, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Number of carbons is 34, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Number of carbons is 9, not 15 '
                                               '(not a sesquiterpene)'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Number of carbons is 26, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Number of carbons is 22, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Number of carbons is 14, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Number of carbons is 32, not '
                                               '15 (not a sesquiterpene)'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Number of carbons is 26, not '
                                               '15 (not a sesquiterpene)'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12C[C@@H](C)\\C=C/C(=O)[C@@](C)(O)[C@H](OC(C)=O)[C@@H](OC(=O)CC(C)C)[C@@]1([H])C(=C)C(=O)O2',
                                      'name': 'neurolenin B',
                                      'reason': 'Number of carbons is 22, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'O1[C@]2([C@@H](O)C[C@]1(O)C(C[C@]3(OC(=O)C([C@]3([C@H](OC(=O)/C(/C)=C\\C)C2)[H])=C)[H])CO)C',
                                      'name': '4,5-Dihydroniveusin A',
                                      'reason': 'Number of carbons is 20, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'O1[C@]2([C@]([C@@H](C1=O)C)([C@@H](OC(=O)C)CC(=CCCC(=C2)C)C)[H])[H]',
                                      'name': 'Acetylbalchanolide',
                                      'reason': 'Number of carbons is 17, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'COC(=O)C1=C/[C@H]2O[C@H]2\\C(C)=C\\[C@H]2OC(=O)C(=C)[C@@H]2[C@H](OC(=O)[C@](C)(O)[C@H](C)OC(C)=O)[C@H]\\1OC(C)=O',
                                      'name': 'Melampodinin',
                                      'reason': 'Number of carbons is 25, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'O1[C@]2([C@]([C@H](C1=O)C)(CCC(=CC[C@H](O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)C(=C2)C)C)[H])[H]',
                                      'name': '(3R,3aS,6E,9S,10E,11aS)-3,6,10-trimethyl-9-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3a,4,5,8,9,11a-hexahydro-3H-cyclodeca[b]furan-2-one',
                                      'reason': 'Number of carbons is 21, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'O1[C@]2([C@@]([C@H](OC)[C@@H](OC(=O)C(C)=C)C(=CCCC(=C2)CO)C=O)(C(C1=O)=C)[H])[H]',
                                      'name': '[(3aS,4S,5S,6E,11aR)-6-formyl-10-(hydroxymethyl)-4-methoxy-3-methylidene-2-oxo-3a,4,5,8,9,11a-hexahydrocyclodeca[b]furan-5-yl] '
                                              '2-methylprop-2-enoate',
                                      'reason': 'Number of carbons is 20, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'O1C2C(C(OC(=O)C(CC)C)C(OC)C(=CCCC(=C2)C)C=O)C(C1=O)=C',
                                      'name': '[(6Z,10E)-6-Formyl-5-methoxy-10-methyl-3-methylidene-2-oxo-3a,4,5,8,9,11a-hexahydrocyclodeca[b]furan-4-yl] '
                                              '2-methylbutanoate',
                                      'reason': 'Number of carbons is 21, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'C[C@H]1O[C@@]1(C)C(=O)O[C@@H]1\\C=C(C)/[C@@H](O)C[C@H](OC(C)=O)\\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12',
                                      'name': 'Eupacunoxin',
                                      'reason': 'Number of carbons is 22, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'CC(=O)O[C@H]1C\\C=C(C)\\C[C@@H](OC(=O)C(\\CO)=C\\CO)[C@@H]2[C@H](OC(=O)C2=C)\\C=C1\\C',
                                      'name': 'Hiyodorilactone A',
                                      'reason': 'Number of carbons is 22, not '
                                                '15 (not a sesquiterpene)'},
                                  {   'smiles': 'CC(=C)C(=O)O[C@@H]1C[C@@]2(C)O[C@@H]2C[C@H](OC(=O)C(C)=C)\\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12',
                                      'name': 'Erioflorin methacrylate',
                                      'reason': 'Number of carbons is 23, not '
                                                '15 (not a sesquiterpene)'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a **germacranolide**, we need '
                 'to check for the following characteristics:\n'
                 '\n'
                 '- **Sesquiterpene**: The molecule should have 15 carbon '
                 'atoms.\n'
                 '- **Germacrane skeleton**: Presence of a 10-membered ring.\n'
                 '- **Lactone group**: A cyclic ester functional group within '
                 'a ring.\n'
                 '\n'
                 'Here is the Python code that implements these checks using '
                 'RDKit:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 24,
    'num_false_positives': 24,
    'num_true_negatives': 142191,
    'num_false_negatives': 61,
    'num_negatives': None,
    'precision': 0.5,
    'recall': 0.2823529411764706,
    'f1': 0.3609022556390977,
    'accuracy': 0.9994026704146171,
    'negative_predictive_value': 0.9995711835334476}