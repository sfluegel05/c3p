"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:37671 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative where one or more oxygens or hydroxyl groups
    are replaced by sulfur or -SR groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a sugar-like structure (ring with multiple hydroxyl groups)
    sugar_pattern = Chem.MolFromSmarts("[C;H1,H2][OH]")  # Pattern for a carbon with a hydroxyl group
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 3:  # At least 3 hydroxyl groups for a sugar-like structure
        return False, "Not enough hydroxyl groups for a sugar-like structure"

    # Check for sulfur atoms in the molecule
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if sulfur is attached to any part of the molecule, not just the sugar backbone
    sulfur_attached = False
    for sulfur in sulfur_atoms:
        neighbors = sulfur.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                sulfur_attached = True
                break
        if sulfur_attached:
            break

    if not sulfur_attached:
        return False, "Sulfur not attached to the molecule"

    # Check for the presence of a ring structure (carbohydrate-like)
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring structure found"

    return True, "Contains a sugar-like structure with sulfur attached"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73754',
                          'name': 'thiosugar',
                          'definition': 'A carbohydrate derivative in which '
                                        'one or more of the oxygens or hydroxy '
                                        'groups of the parent carbohydrate is '
                                        'replaced by sulfur or -SR, where R '
                                        'can be hydrogen or any group.',
                          'parents': ['CHEBI:33261', 'CHEBI:63299'],
                          'xrefs': ['PMID:16240117', 'PMID:23330717'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: SMILES: '
               'S([C@@H](O)C[C@H]1[C@H](O[C@H]2O[C@H](C[C@H]([C@H]2O)N(C)C)C)[C@@H]([C@@H](O)CC(=O)O[C@H]([C@H](C)[C@H]3[C@@](C=CC([C@@H](C1)C)=O)(O3)C)CC)C)C[C@H](N)C(=O)O '
               'NAME: 6108C REASON: WRONGLY CLASSIFIED Contains a sugar-like '
               'structure with sulfur attached\n'
               ' * SMILES: '
               'CC(=O)N[C@@H](CSCO)C(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O '
               'NAME: S-(hydroxymethyl)mycothiol REASON: WRONGLY CLASSIFIED '
               'Contains a sugar-like structure with sulfur attached\n'
               ' * SMILES: '
               'CSC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxy-3-(methylthio)propanoyl-CoA REASON: WRONGLY '
               'CLASSIFIED Contains a sugar-like structure with sulfur '
               'attached\n'
               'False negatives: SMILES: '
               'S([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)C(=NO)CCCCCCSC '
               'NAME: 6-(methylthio)hexyldesulfoglucosinolate REASON: MISSED '
               'Sulfur not attached to the sugar backbone\n'
               ' * SMILES: '
               '[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)SC(=NO)CC=2C=CC=CC2 '
               'NAME: desulfoglucotropeolin REASON: MISSED Sulfur not attached '
               'to the sugar backbone\n'
               ' * SMILES: '
               '[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/C(CC)C '
               'NAME: glucocochlearin REASON: MISSED Sulfur not attached to '
               'the sugar backbone\n'
               ' * SMILES: '
               'CC[P+](CC)(CC)[Au-]S[C@@H]1O[C@H](COC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@H]1OC(C)=O '
               'NAME: auranofin REASON: MISSED Not enough hydroxyl groups for '
               'a sugar-like structure\n'
               ' * SMILES: '
               '[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCC=C '
               'NAME: gluconapin REASON: MISSED Sulfur not attached to the '
               'sugar backbone\n'
               ' * SMILES: '
               'CCOC1=CC=C(CC2=CC(=CC=C2Cl)[C@@H]2O[C@H](SC)[C@@H](O)[C@H](O)[C@H]2O)C=C1 '
               'NAME: sotagliflozin REASON: MISSED Sulfur not attached to the '
               'sugar backbone\n'
               ' * SMILES: '
               'O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1SC(CCCCCCCS(C)=O)=NOS(O)(=O)=O)O)O)O)CO '
               'NAME: glucoibarin REASON: MISSED Sulfur not attached to the '
               'sugar backbone\n'
               ' * SMILES: '
               'CCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1S[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O '
               'NAME: propyl '
               'beta-D-mannopyranosyl-(1->2)-beta-D-mannopyranosyl-(1->2)-2-thio-beta-D-mannopyranoside '
               'REASON: MISSED Sulfur not attached to the sugar backbone\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](SC(Cc2c[nH]c3cccc(O)c23)=NOS(O)(=O)=O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: 4-hydroxyglucobrassicin REASON: MISSED Sulfur not '
               'attached to the sugar backbone\n'
               ' * SMILES: '
               'COn1cc(CC(S[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)=NOS(O)(=O)=O)c2ccccc12 '
               'NAME: neoglucobrassicin REASON: MISSED Sulfur not attached to '
               'the sugar backbone\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1=CC=C2C(=C1)C3=C(C4=C(C=CC=N4)C=C3)NS2(=O)=O',
                                     'name': '5H-quinolino[8,7-c][1,2]benzothiazine '
                                             '6,6-dioxide',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'O([C@@H]1C([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@]([C@H](O)C4)(CCC(C5)(C)C)C(O)=O)[H])=CC3)C)(CC2)C)[H])(CC1)C)[H])(C)C)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO',
                                     'name': 'Echinocystic acid 3-glucoside',
                                     'reason': 'No sulfur atoms found'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=C(C=C2)NS(=O)(=O)C3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C',
                                     'name': 'N-[(4R,7R,8S)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]benzenesulfonamide',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'O=C1O[C@](C(=O)OC)(CC2=CC(=C(O)C=C2)CC(=O)C)C(=C1O)C3=CC=C(OC)C=C3',
                                     'name': 'Versicolactone A',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'O=C(N[C@@H](CO)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1NC=NC1)CCC(O)=O',
                                     'name': 'His-Glu-Ser',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'O=C1N(O)CCOCCNC(=O)CCC(=O)N(O)CCOCCNC(CCC(N(CCOCCNC(CC1)=O)O)=O)=O',
                                     'name': 'Desferrioxamine Et3',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O',
                                     'name': 'PE(18:3(9Z,12Z,15Z)/18:1(11Z))',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'S(C[C@H](N)C(=O)N[C@H]1[C@@H](O)C(N2C3=NC=NC(=C3N=C2)N(C)C)OC1CO)C',
                                     'name': 'Cystocin',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'COc1cc(O)c(CC=C(C)C)cc1C(=O)[C@@H](O)Cc1ccc(O)cc1',
                                     'name': 'lespeflorin C3',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'},
                                 {   'smiles': 'NC(=O)OP(O)(O)=O',
                                     'name': 'carbamoyl phosphate',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar-like structure'}],
    'sample_false_negatives': [   {   'smiles': 'CC[P+](CC)(CC)[Au-]S[C@@H]1O[C@H](COC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@H]1OC(C)=O',
                                      'name': 'auranofin',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': 'CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O',
                                      'name': '1-(methylthio)ribulose '
                                              '5-phosphate',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': 'ClC(C(NC(=O)C1N(CC(C1)CCC)C)C2OC(SC)C(OP(=O)(O)O)C(C2O)O)C',
                                      'name': 'Clindamycin 2-Phosphate',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': 'S(C1OC(C(OC(=O)C)C(OC(=O)C)C1N2C(=O)C=3C(C2=O)=CC=CC3)COC(=O)C)C',
                                      'name': 'Methyl '
                                              '3,4,6-tri-O-acetyl-2-deoxy-2-(1,3-dioxo-1,3-dihydro-2H-isoindol-2-yl)-1-thiohexopyranoside',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': 'S([C@H]1O[C@H](CN2C=3N=CNC(C3N=C2)=O)[C@H]([C@H]1O)O)C',
                                      'name': '9-[[(2R,3S,4R,5R)-3,4-dihydroxy-5-methylsulanyloxolan-2-yl]methyl]-1H-purin-6-one',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': '[C@](NC(=O)[C@@H]1C[C@H](CN1C)CCC)([C@H](C)Cl)([C@]2(O[C@H](SC)[C@@H]([C@H]([C@H]2O)O)OC(CCCCCCCCCCCCCCC)=O)[H])[H]',
                                      'name': 'clindamycin palmitate',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'},
                                  {   'smiles': '[C@@]1([C@@H]([C@@H]([C@H]([C@H](O1)SC)OC(CCCCCCCCCCCCCCC)=O)O)O)([C@](NC([C@H]2N(C[C@@H](C2)CCC)C)=O)([C@@H](Cl)C)[H])[H].Cl',
                                      'name': 'Clindamycin palmitate '
                                              'hydrochloride',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar-like structure'}],
    'attempt': 1,
    'reasoning': 'The previous code had several issues that led to incorrect '
                 'classifications:\n'
                 '\n'
                 '1. **False Positives**: The code incorrectly classified '
                 'molecules that had sulfur attached to a sugar-like structure '
                 'but were not thiosugars. This happened because the code only '
                 'checked for the presence of sulfur and hydroxyl groups '
                 'without considering the specific structural requirements of '
                 'thiosugars.\n'
                 '\n'
                 '2. **False Negatives**: The code missed several thiosugars '
                 'because it required sulfur to be directly attached to the '
                 'sugar backbone. However, thiosugars can have sulfur attached '
                 'to other parts of the molecule, such as side chains or '
                 'substituents.\n'
                 '\n'
                 '3. **Insufficient Structural Checks**: The code did not '
                 'adequately check for the presence of a carbohydrate-like '
                 'structure. It only checked for hydroxyl groups but did not '
                 'verify the presence of a ring structure or the typical '
                 'arrangement of hydroxyl groups in a carbohydrate.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '- Ensure that the molecule has a carbohydrate-like '
                 'structure, including a ring and multiple hydroxyl groups.\n'
                 '- Allow sulfur to be attached to any part of the molecule, '
                 'not just the sugar backbone.\n'
                 '- Check for the presence of sulfur or -SR groups in the '
                 'molecule.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 61,
    'num_false_positives': 100,
    'num_true_negatives': 13894,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.37888198757763975,
    'recall': 0.8970588235294118,
    'f1': 0.5327510917030568,
    'accuracy': 0.99239084056322,
    'negative_predictive_value': 0.9994964391051003}