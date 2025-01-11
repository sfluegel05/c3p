"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:37671 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxamidine pattern: RC(=NR)NR2
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3H2,NX3H1,NX3H0]")
    if not mol.HasSubstructMatch(carboxamidine_pattern):
        return False, "No carboxamidine structure (RC(=NR)NR2) found"

    return True, "Contains the carboxamidine structure (RC(=NR)NR2)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35359',
                          'name': 'carboxamidine',
                          'definition': 'Compounds having the structure '
                                        'RC(=NR)NR2. The term is used as a '
                                        'suffix in systematic nomenclature to '
                                        'denote the -C(=NH)NH2 group including '
                                        'its carbon atom.',
                          'parents': ['CHEBI:2634', 'CHEBI:35352'],
                          'xrefs': ['KEGG:C06060'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': '[Si](OC(=O)CCCCCCCCCCC)(C)(C)C',
                                     'name': 'Dodecanoic acid, trimethylsilyl '
                                             'ester',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)NC4=C(C=C(C=C4)OC)OC)[C@H](C)CO',
                                     'name': 'LSM-30778',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'C(C([C@H]([C@@H]([C@H]([C@@H](CO)O)O)O)O)=O)O',
                                     'name': 'D-ido-heptulose',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'O=C1C=C(CC(O)C)O[C@]1(/C=C/C=C/CC)C',
                                     'name': 'Terrefuranone',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'C1CC1CNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@H](O2)CO)OC4=C3C=C(C=C4)NC(=O)CC5=CC=NC=C5',
                                     'name': '2-[(1R,3S,4aS,9aR)-1-(hydroxymethyl)-6-[(1-oxo-2-pyridin-4-ylethyl)amino]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-(cyclopropylmethyl)acetamide',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]4CO)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O[C@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H]9O)CO)[C@H](O)[C@H]8NC(=O)C)CO)[C@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)[C@H](O)[C@@H]%11O[C@@H]%14O[C@@H]([C@@H](O[C@@H]%15O[C@@H]([C@H](O)[C@H](O)[C@H]%15O)CO)[C@H](O[C@@H]%16O[C@H]([C@@H](O)[C@@H](O)[C@@H]%16O)C)[C@H]%14NC(=O)C)CO)CO',
                                     'name': 'CID 91853359',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)C2=NN3C(=C2)N=C4CCCC4=C3N5CCC6(CC5)OCCO6',
                                     'name': 'LSM-28273',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': '[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O',
                                     'name': '4-deoxy-4-sulfo-D-erythrose(1-)',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': 'O[C@@H](CCCCC)\\C=C\\CCCCCCCCCC(O)=O',
                                     'name': '13S-HOME(11E)',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'},
                                 {   'smiles': '[H][C@]1(O[C@@](C[C@H](O)[C@H]1NC(=O)CO)(O[C@H]1[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]2[C@@H](NC(C)=O)[C@H](O)O[C@H](CO)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@@H]1O)C(O)=O)[C@H](O)[C@H](O)CO',
                                     'name': 'alpha-Neup5Gc-(2->3)-beta-D-Galp-(1->3)-[alpha-L-Fucp-(1->4)]-beta-D-GlcNAcp',
                                     'reason': 'No carboxamidine structure '
                                               '(RC(=NR)NR2) found'}],
    'sample_false_negatives': [   {   'smiles': 'CNC(NCCSCc1csc(CN(C)C)n1)=C[N+]([O-])=O',
                                      'name': 'nizatidine',
                                      'reason': 'No carboxamidine structure '
                                                '(RC(=NR)NR2) found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as carboxamidines, we need to '
                 'identify the presence of the `RC(=NR)NR2` structure, where R '
                 'can be any substituent. The key feature is the presence of a '
                 'carbon double-bonded to a nitrogen (C=NR) and single-bonded '
                 'to another nitrogen (C-NR2). This can be detected using '
                 'SMARTS patterns in RDKit.\n'
                 '\n'
                 'Here’s the Python program to classify carboxamidines:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 38,
    'num_false_positives': 100,
    'num_true_negatives': 4740,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.2753623188405797,
    'recall': 0.9743589743589743,
    'f1': 0.4293785310734463,
    'accuracy': 0.9792990366878459,
    'negative_predictive_value': 0.9997890740350137}