"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    These are steroids with a ketone at position 3 and a double bond between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid ring system (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) that's part of the first ring of the steroid
    ketone_pattern = Chem.MolFromSmarts("[#6]~1~[#6]=,:[#6]~[#6](=[O])~[#6]~[#6]~[#6]~2")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone found at position 3"

    # Check for double bond between positions 1 and 2
    # This pattern specifically looks for the double bond between C1 and C2,
    # ensuring it's part of the ring with the ketone at C3
    double_bond_pattern = Chem.MolFromSmarts("[#6]~1~[#6]=,:[#6]~[#6](=[O])~[#6]~[#6]~[#6]~2")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond between positions 1 and 2"

    # Additional check to ensure the molecule has reasonable size for a steroid
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Steroids typically have more than 20 atoms
        return False, "Molecule too small to be a steroid"

    return True, "Molecule contains steroid core with 3-oxo group and Delta(1) double bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:20156',
                          'name': '3-oxo-Delta(1) steroid',
                          'definition': 'Any 3-oxo steroid that contains a '
                                        'double bond between positions 1 and '
                                        '2.',
                          'parents': ['CHEBI:47788', 'CHEBI:51689'],
                          'xrefs': [   'KEGG:C02941',
                                       'MetaCyc:3-Oxo-Delta-1-Steroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'No steroid core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)COC(=O)CCC(O)=O)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'prednisolone succinate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C=CC2=CC(=O)C=C[C@]12C',
                                      'name': 'Delta(6)-prednisone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@H]([C@H](O)[C@@H]1OC(=O)[C@@H](C)[C@@H]1C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': 'paraminabeolide E',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CC[C@@]3([H])[C@@]4(C)C=CC(=O)[C@@H](C)[C@]4([H])[C@H](OC(C)=O)C(=O)[C@]3(C)[C@@]1(C)C[C@H](OC(C)=O)\\C2=C(\\CCCC(C)(C)O)C(O)=O',
                                      'name': '6beta,16beta-diacetoxy-25-hydroxy-3,7-dioxy-29-nordammara-1,17(20)-dien-21-oic '
                                              'acid',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@H](C)[C@](OC(=O)CCCC)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'betamethasone valerate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)([C@H](O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': '20alpha-dihydroprednisone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)C=C[C@]4(C)[C@@]3([H])CC[C@]12COC(C)=O)[C@H](C)[C@H]1CC(C)=C(C)C(=O)O1',
                                      'name': 'minabeolide 2',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)CC=C1[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'vamorolone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@H](C)[C@](OC(=O)C(C)C)(C(=O)COC(C)=O)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'betamethasone acibutate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@H](C)[C@](OC(=O)CC)(C(=O)CCl)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'clobetasol propionate',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify 3-oxo-Delta(1) '
                 'steroids. The key characteristics are:\n'
                 '1. Must be a steroid (four fused rings)\n'
                 '2. Must have a ketone (=O) at position 3\n'
                 '3. Must have a double bond between positions 1 and 2\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142226,
    'num_false_negatives': 74,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9994799718903724,
    'negative_predictive_value': 0.9994799718903724}