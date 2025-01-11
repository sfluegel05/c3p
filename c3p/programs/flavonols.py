"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:28413 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are hydroxyflavones with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonoid core with 3-OH group and 4-ketone (flavonol backbone)
    # [#6]1=[#6]-c2c([#6](=[O])[#6]([OH])=[#6]1)c([OH])cc([OH])c2
    flavonol_pattern = Chem.MolFromSmarts('[#6]1=[#6]-c2c([#6](=[O])[#6]([OH])=[#6]1)c([OH,O-1])cc([OH,O-1])c2')
    
    # Alternative pattern to catch variations
    flavonol_pattern2 = Chem.MolFromSmarts('O=C1C(O)=C(Oc2ccccc12)c1ccccc1')
    
    if not (mol.HasSubstructMatch(flavonol_pattern) or mol.HasSubstructMatch(flavonol_pattern2)):
        return False, "Missing flavonol core structure (3-hydroxyflavone backbone)"

    # Check for ketone at position 4
    ketone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2OC(c2ccccc2)=C1O')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group at position 4"

    # Count oxygen atoms (should have at least 4: 3-OH, 4=O, and ring O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Insufficient oxygen atoms for flavonol structure"

    # Check for aromatic rings (should have at least 2)
    aromatic_rings = 0
    for atom in mol.GetAtoms():
        if atom.IsInRing() and atom.GetIsAromatic():
            aromatic_rings += 1
    if aromatic_rings < 8:  # Each aromatic ring contributes multiple aromatic atoms
        return False, "Missing required aromatic rings"

    # Additional check for 2,3-double bond
    double_bond_pattern = Chem.MolFromSmarts('C1=COc2ccccc2C1=O')
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "Missing required double bond between positions 2 and 3"

    return True, "Contains 3-hydroxyflavone backbone with required structural features"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28802',
                          'name': 'flavonols',
                          'definition': 'Any hydroxyflavone in which is the '
                                        'ring hydrogen at position 3 of the '
                                        'heterocyclic ring is replaced by a '
                                        'hydroxy group.',
                          'parents': ['CHEBI:192499', 'CHEBI:24698'],
                          'xrefs': ['MetaCyc:Flavonols', 'Wikipedia:Flavonol'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Missing flavonol core '
                                               'structure (3-hydroxyflavone '
                                               'backbone)'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)=CCC\\C(C)=C\\Cc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1',
                                      'name': '8-geranylchrysin',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'COc1cc(OC)c2c(c1)oc(-c1ccc(OC)c(OC)c1)c(O)c2=O',
                                      'name': "quercetin 5,7,3',4'-tetramethyl "
                                              'ether',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'COc1ccc(cc1)-c1oc2cc(O)cc(O)c2c(=O)c1O',
                                      'name': 'kaempferide',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)ccc2c1oc(-c1ccc(O)c(O)c1)c(O)c2=O',
                                      'name': 'fisetin 8-C-glucoside',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'CC(C)=CCc1c(O)c(CC=C(C)C)c2oc(cc(=O)c2c1O)-c1ccccc1',
                                      'name': '6,8-di-(3,3-dimethylallyl)chrysin',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'CC(C)=CCc1c(O)cc2oc(-c3ccccc3)c(O)c(=O)c2c1O',
                                      'name': '6-(3,3-dimethylallyl)galangin',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2ccc(cc2)-c2oc3cc(O)cc(O)c3c(=O)c2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': 'kaempferol '
                                              "4'-O-beta-D-glucopyranoside",
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'COc1cc(ccc1O)-c1oc2cc(O)cc(O)c2c(=O)c1O',
                                      'name': 'isorhamnetin',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'Oc1cc(O)c2c(c1)oc(-c1cc(O)c(O)c(O)c1)c(O)c2=O',
                                      'name': 'myricetin',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc(O)c3c(c2)oc(-c2ccc(O)cc2)c(O)c3=O)[C@H](O)[C@@H](O)[C@H]1O',
                                      'name': 'kaempferol '
                                              '7-O-beta-D-galactopyranoside',
                                      'reason': 'Missing flavonol core '
                                                'structure (3-hydroxyflavone '
                                                'backbone)'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify flavonols. '
                 'Flavonols are a subclass of flavonoids characterized by '
                 'having a 3-hydroxyflavone backbone. The key structural '
                 'features are:\n'
                 '\n'
                 '1. A C6-C3-C6 skeleton (basic flavonoid structure)\n'
                 '2. A ketone group at position 4\n'
                 '3. A hydroxyl group at position 3\n'
                 '4. A double bond between positions 2 and 3\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142236,
    'num_false_negatives': 64,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.999550245959241,
    'negative_predictive_value': 0.999550245959241}