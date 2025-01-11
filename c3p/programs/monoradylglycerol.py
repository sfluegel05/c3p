"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
A lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic property checks
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight < 100 or mol_weight > 500:
        return False, f"Molecular weight {mol_weight:.1f} outside expected range (100-500)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:
        return False, f"Too few carbons ({c_count}) for monoradylglycerol"
    if not (3 <= o_count <= 4):
        return False, f"Invalid number of oxygen atoms ({o_count}), expected 3-4"

    # Look for glycerol backbone - more general pattern
    glycerol_pattern = Chem.MolFromSmarts("[OX2,OH1][CH2][CH]([OX2,OH1])[CH2][OX2,OH1]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count free hydroxyls
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches != 2:
        return False, f"Found {oh_matches} free hydroxyl groups, need exactly 2"

    # Define substituent patterns
    patterns = {
        'acyl': [
            # 1/3-acyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CX3](=[OX1])[#6]",
            # 2-acyl
            "[OH1][CH2][CH]([OX2][CX3](=[OX1])[#6])[CH2][OH1]"
        ],
        'alkyl': [
            # 1/3-alkyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CX4;!$(C=O)]",
            # 2-alkyl
            "[OH1][CH2][CH]([OX2][CX4;!$(C=O)])[CH2][OH1]"
        ],
        'alkenyl': [
            # 1/3-alkenyl
            "[OX2]([CH2][CH]([OH1])[CH2][OH1])[CH2][CH]=[CH]",
            # 2-alkenyl
            "[OH1][CH2][CH]([OX2][CH2][CH]=[CH])[CH2][OH1]"
        ]
    }

    # Count substituents
    total_matches = 0
    substituent_type = None
    
    for subst_type, subst_patterns in patterns.items():
        for pattern in subst_patterns:
            matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
            if matches > 0:
                total_matches += matches
                substituent_type = subst_type

    if total_matches != 1:
        return False, f"Found {total_matches} substituents, need exactly 1"

    # Additional validation for chain length
    chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    chain_matches = len(mol.GetSubstructMatches(chain_pattern))
    if chain_matches < 1:
        return False, "Substituent chain too short"

    return True, f"Valid monoradylglycerol with one {substituent_type} substituent"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76575',
                          'name': 'monoradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing a '
                                        'single acyl, alkyl or alk-1-enyl '
                                        'substituent at an unspecified '
                                        'position.',
                          'parents': ['CHEBI:35741'],
                          'xrefs': ['LIPID_MAPS_class:LMGL0101'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6373626373626374 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OCC(O)CO '
               'NAME: 1-monolinolein REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: O=C(OC[C@@H](O)CO)CCCCCCCCCCCCCCCC(C)C NAME: '
               'AKD-2B1 REASON: CORRECT Valid monoradylglycerol with one acyl '
               'substituent\n'
               ' * SMILES: '
               'O(C(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO '
               'NAME: MG(22:4(7Z,10Z,13Z,16Z)/0:0/0:0) REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCC(=O)OCC(O)CO NAME: '
               'rac-1-monodecanoylglycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: OCC(COC(=O)CCCCCCC/C=C\\CCCCCC)O NAME: '
               '1-[(9Z)-hexadecenoyl]glycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCCC)C[C@@H](O)CO NAME: '
               'MG(15:0/0:0/0:0) REASON: CORRECT Valid monoradylglycerol with '
               'one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-myristoyl-sn-glycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-decanoyl-sn-glycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: '
               'O(C[C@@H](O)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC '
               'NAME: MG(22:5(4Z,7Z,10Z,13Z,16Z)/0:0/0:0) REASON: CORRECT '
               'Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-palmitoyl-sn-glycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCC/C=C\\CCCCCCCC)C[C@@H](O)CO '
               'NAME: MG(24:1(15Z)/0:0/0:0) REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: C(CCCCCCCCCCCC)CCCCCCCCC(OCC(O)CO)=O NAME: '
               '1-behenoylglycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C/CCCCCC)CC(O)CO NAME: '
               'MG(18:1(11E)/0:0/0:0)[rac] REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O(C(=O)CCCCCCCCCC=CCCCCCCCC)CC(O)CO NAME: '
               '2,3-dihydroxypropylicos-11-enoate REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(=O)OCC(O)CO NAME: '
               'rac-1-monopalmitoylglycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCC(=O)OCC(O)CO NAME: 1-monooctanoylglycerol '
               'REASON: CORRECT Valid monoradylglycerol with one acyl '
               'substituent\n'
               ' * SMILES: '
               'C(\\C/C=C\\C/C=C\\C/C=C\\CC)=C\\C/C=C\\C/C=C\\CCC(=O)OCC(CO)O '
               'NAME: 1-(all-cis-docosa-4,7,10,13,16,19-hexaenoyl)-glycerol '
               'REASON: CORRECT Valid monoradylglycerol with one acyl '
               'substituent\n'
               ' * SMILES: OCC(COC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)O NAME: '
               '1-nervonoylglycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO NAME: '
               '1-stearoyl-sn-glycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C NAME: '
               '2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate REASON: '
               'CORRECT Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-stearoyl-sn-glycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: CCCCCCCC(=O)OC[C@H](O)CO NAME: '
               '3-octanoyl-sn-glycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(O)CO NAME: '
               '1-oleoylglycerol REASON: CORRECT Valid monoradylglycerol with '
               'one acyl substituent\n'
               ' * SMILES: CCCCCCCCCC(=O)OC[C@@H](O)CO NAME: '
               '1-decanoyl-sn-glycerol REASON: CORRECT Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OCC(O)CO NAME: '
               'rac-1-monomyristoylglycerol REASON: CORRECT Valid '
               'monoradylglycerol with one acyl substituent\n'
               'False positives: SMILES: C(CCCCCCCC)CCCOCC(O)CO NAME: '
               '1-O-dodecylglycerol REASON: WRONGLY CLASSIFIED Valid '
               'monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: '
               'O(C[C@H](OC)C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO '
               'NAME: '
               '1-o-(2r-methoxy-4z,7z,10z,13z,16z-docosapentaenyl)-sn-glycerol '
               'REASON: WRONGLY CLASSIFIED Valid monoradylglycerol with one '
               'alkyl substituent\n'
               ' * SMILES: O(CCCCCCCCCCCCCC)CC(O)CO NAME: '
               '1-O-Tetradecylglycerol REASON: WRONGLY CLASSIFIED Valid '
               'monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: O(C1C(CCC(C1)C)C(C)C)CC(O)CO NAME: '
               '3-[[5-Methyl-2-(1-methylethyl)cyclohexyl]oxy]-1,2-propanediol '
               'REASON: WRONGLY CLASSIFIED Valid monoradylglycerol with one '
               'alkyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCOCC(O)CO NAME: batilol REASON: '
               'WRONGLY CLASSIFIED Valid monoradylglycerol with one alkyl '
               'substituent\n'
               ' * SMILES: O(C[C@H](OC)C/C=C\\CCCCCCCCCC)C[C@@H](O)CO NAME: '
               '1-O-(2R-methoxy-4Z-pentadecenyl)-sn-glycerol REASON: WRONGLY '
               'CLASSIFIED Valid monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: OCC(O)COC(=O)c1ccccc1Nc1ccnc2cc(Cl)ccc12 NAME: '
               'glafenine REASON: WRONGLY CLASSIFIED Valid monoradylglycerol '
               'with one acyl substituent\n'
               ' * SMILES: '
               'O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C[C@@H](O)CO NAME: '
               'MG(20:4(8Z,11Z,14Z,17Z)/0:0/0:0) REASON: WRONGLY CLASSIFIED '
               'Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: '
               'C(F)(F)(F)C1=C2C(C(NC3=C(C(OCC(CO)O)=O)C=CC=C3)=CC=N2)=CC=C1 '
               'NAME: Floctafenine REASON: WRONGLY CLASSIFIED Valid '
               'monoradylglycerol with one acyl substituent\n'
               ' * SMILES: '
               'O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C[C@@H](O)CO NAME: '
               'MG(18:4(6Z,9Z,12Z,15Z)/0:0/0:0) REASON: WRONGLY CLASSIFIED '
               'Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O(CCC\\C=C\\C(=O)CCCCCCCCCC)C[C@@H](O)CO NAME: '
               'Ceratodictyol A REASON: WRONGLY CLASSIFIED Valid '
               'monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO '
               'NAME: MG(20:2(11Z,14Z)/0:0/0:0) REASON: WRONGLY CLASSIFIED '
               'Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O(CC(O)CO)C NAME: 3-Methoxy-1,2-propanediol REASON: '
               'WRONGLY CLASSIFIED Valid monoradylglycerol with one alkyl '
               'substituent\n'
               ' * SMILES: O([C@H](CCCCCCCCCCCCCC)COC[C@@H](O)CO)C NAME: '
               '1-O-(2R-methoxy-hexadecyl)-sn-glycerol REASON: WRONGLY '
               'CLASSIFIED Valid monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: OCC(O)COC(=O)CS NAME: glyceryl monothioglycolate '
               'REASON: WRONGLY CLASSIFIED Valid monoradylglycerol with one '
               'acyl substituent\n'
               ' * SMILES: O([C@H](C/C=C\\CCCCCCCCCCCCCCC)COC[C@@H](O)CO)C '
               'NAME: 1-O-(2R-methoxy-4Z-eicosenyl)-sn-glycerol REASON: '
               'WRONGLY CLASSIFIED Valid monoradylglycerol with one alkyl '
               'substituent\n'
               ' * SMILES: O(CCCCCC)C[C@@H](O)CO NAME: 1-Hexylglycerol REASON: '
               'WRONGLY CLASSIFIED Valid monoradylglycerol with one alkyl '
               'substituent\n'
               ' * SMILES: O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO '
               'NAME: MG(20:3(8Z,11Z,14Z)/0:0/0:0) REASON: WRONGLY CLASSIFIED '
               'Valid monoradylglycerol with one acyl substituent\n'
               ' * SMILES: O([C@H](C/C=C\\CCCCCCCCCCCC)COC[C@@H](O)CO)C NAME: '
               '1-O-(2R-methoxy-4Z-heptadecenyl)-sn-glycerol REASON: WRONGLY '
               'CLASSIFIED Valid monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: O(C(CCCCCCCCCCCCCC)COC[C@@H](O)CO)C NAME: '
               '1-O-(2-methoxy-hexadecyl)-sn-glycerol REASON: WRONGLY '
               'CLASSIFIED Valid monoradylglycerol with one alkyl substituent\n'
               ' * SMILES: CCCCCCCCCCCCCCCCOCC(O)CO NAME: 1-O-palmitylglycerol '
               'REASON: WRONGLY CLASSIFIED Valid monoradylglycerol with one '
               'alkyl substituent\n'
               ' * SMILES: CCCCCCCC\\C=C/CCCCCCCCOCC(O)CO NAME: '
               '(9Z)-1-O-octadec-9-enyl glycerol REASON: WRONGLY CLASSIFIED '
               'Valid monoradylglycerol with one alkyl substituent\n'
               'False negatives: SMILES: O(C(=O)CCCCCCCC=CCC=CCCCCC)C(CO)CO '
               'NAME: 2-Linoleoyl glycerol REASON: MISSED Found 0 '
               'substituents, need exactly 1\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/15:0/0:0) REASON: MISSED Found 0 substituents, need '
               'exactly 1\n'
               ' * SMILES: CCCCCCCCCC(=O)OC(CO)CO NAME: 2-decanoylglycerol '
               'REASON: MISSED Found 0 substituents, need exactly 1\n'
               ' * SMILES: '
               'O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC NAME: '
               'MG(0:0/22:5(4Z,7Z,10Z,13Z,16Z)/0:0) REASON: MISSED Found 0 '
               'substituents, need exactly 1\n'
               ' * SMILES: O=C1OC(/C(=C/C(=O)OC[C@H](O)CO)/C)=CC(=C1)OC NAME: '
               'Pestalotiopyrone I REASON: MISSED Invalid number of oxygen '
               'atoms (7), expected 3-4\n'
               ' * SMILES: O(C(CO)CO)C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: MG(0:0/18:4(6Z,9Z,12Z,15Z)/0:0) REASON: MISSED Found 0 '
               'substituents, need exactly 1\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC(CO)CO NAME: '
               '2-linoleoylglycerol REASON: MISSED Found 0 substituents, need '
               'exactly 1\n'
               ' * SMILES: O1C(C1C/C=C\\CCCCC)CCCCCCCC(OC[C@@H](O)CO)=O NAME: '
               'MG(18:1(12Z)-O(9S,10R)/0:0/0:0) REASON: MISSED Invalid number '
               'of oxygen atoms (5), expected 3-4\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/22:0/0:0) REASON: MISSED Found 0 substituents, need '
               'exactly 1\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(=O)OC(CO)CO NAME: '
               '2-palmitoylglycerol REASON: MISSED Found 0 substituents, need '
               'exactly 1\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC(=O)CCNC(O)=O',
                                     'name': 'N-carboxy-beta-alanine',
                                     'reason': 'Too few carbons (4) for '
                                               'monoradylglycerol'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)NC3=CC=CC=C3F)C(=O)N(C[C@H]1OC)C)C)CCC(F)(F)F',
                                     'name': '1-(2-fluorophenyl)-3-[(5S,6R,9R)-5-methoxy-3,6,9-trimethyl-2-oxo-8-(3,3,3-trifluoropropyl)-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'Molecular weight 554.3 outside '
                                               'expected range (100-500)'},
                                 {   'smiles': 'O(C1=CC(=C(O)C=C1)C2=C(O)C=CC(=C2)OC)C',
                                     'name': '5,5′-dimethoxybiphenyl-2,2′-diol',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'COC1=CC=C(CCNC[C@H](O)C2=CC=C(O)C=C2)C=C1OC',
                                     'name': 'denopamine',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': '[Zn++].[S-]C(=S)NCCNC([S-])=S',
                                     'name': 'zineb',
                                     'reason': 'Too few carbons (4) for '
                                               'monoradylglycerol'},
                                 {   'smiles': 'P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCC(C)C)COC(=O)CCCCCCCCC(CC)C)(O)=O)(O)=O',
                                     'name': 'CDP-DG(a-13:0/i-15:0)',
                                     'reason': 'Molecular weight 897.5 outside '
                                               'expected range (100-500)'},
                                 {   'smiles': 'O=C1OC(O)C2=C1C[C@](O)([C@H]3CC(C[C@H]3[C@@H]2O)(C)C)C',
                                     'name': 'Lactarolide A',
                                     'reason': 'Invalid number of oxygen atoms '
                                               '(5), expected 3-4'},
                                 {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'Invalid number of oxygen atoms '
                                               '(5), expected 3-4'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-2,5-Dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-3-yl]acetamide',
                                     'reason': 'Molecular weight 545.2 outside '
                                               'expected range (100-500)'},
                                 {   'smiles': 'COc1ccc(\\C=C/C2CCC=CC2c2ccc(OC)c(OC)c2)cc1OC',
                                     'name': '3-(3,4-dimethoxyphenyl)-4-[(Z)-3,4-dimethoxystyryl]cyclohex-1-ene',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1OC(/C(=C/C(=O)OC[C@H](O)CO)/C)=CC(=C1)OC',
                                      'name': 'Pestalotiopyrone I',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (7), expected 3-4'},
                                  {   'smiles': 'C(\\C/C=C\\C/C=C\\C/C=C\\CC)=C\\C/C=C\\C/C=C\\CCC(=O)OCC(CO)O',
                                      'name': '1-(all-cis-docosa-4,7,10,13,16,19-hexaenoyl)-glycerol',
                                      'reason': 'Substituent chain too short'},
                                  {   'smiles': 'O1C(C1C/C=C\\CCCCC)CCCCCCCC(OC[C@@H](O)CO)=O',
                                      'name': 'MG(18:1(12Z)-O(9S,10R)/0:0/0:0)',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (5), expected 3-4'},
                                  {   'smiles': 'O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO',
                                      'name': '12-HPETE 2-glyceryl ester',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (6), expected 3-4'},
                                  {   'smiles': 'O1C(CCC1CC(=O)C)C(C)C(OCC(O)CO)=O',
                                      'name': '2,3-dihydroxypropyl '
                                              '2-[5-(2-oxopropyl)oxolan-2-yl]propanoate',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (6), expected 3-4'},
                                  {   'smiles': 'O(C[C@@H](O)CO)C(=O)C',
                                      'name': '(R)-glycerol 1-acetate',
                                      'reason': 'Substituent chain too short'},
                                  {   'smiles': 'O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC',
                                      'name': 'MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)',
                                      'reason': 'Substituent chain too short'},
                                  {   'smiles': 'O(CC(O)CO)C(=O)CC',
                                      'name': 'Glycerol 1-propanoate',
                                      'reason': 'Substituent chain too short'},
                                  {   'smiles': '[C@@H]1(/C=C/[C@@H](O)CCCCC)[C@H]([C@H](C[C@H]1O)O)C/C=C\\CCCC(=O)OCC(CO)O',
                                      'name': 'prostaglandin F2alpha '
                                              '1-glyceryl ester',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (7), expected 3-4'},
                                  {   'smiles': 'C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO',
                                      'name': 'prostaglandin D2 2-glyceryl '
                                              'ester',
                                      'reason': 'Invalid number of oxygen '
                                                'atoms (7), expected 3-4'}],
    'attempt': 3,
    'reasoning': 'Let me analyze the issues from the previous attempt:\n'
                 '\n'
                 '1. False Negatives Analysis:\n'
                 '- Many 2-substituted glycerols were missed (e.g., '
                 '"2-Linoleoyl glycerol", "2-decanoylglycerol")\n'
                 '- The glycerol backbone patterns were too specific and '
                 'rigid\n'
                 '- The substituent counting logic was failing to detect valid '
                 'substitutions\n'
                 '\n'
                 '2. False Positives Analysis:\n'
                 '- The program incorrectly classified several alkyl '
                 'glycerols\n'
                 '- Some compounds with complex substituents were wrongly '
                 'identified\n'
                 "- The oxygen count check wasn't sufficient to filter out "
                 'invalid structures\n'
                 '\n'
                 "Here's how we can improve:\n"
                 '\n'
                 '1. Simplify the glycerol backbone detection:\n'
                 '- Use a more general pattern that captures all three '
                 'possible substitution positions\n'
                 '- Focus on the core C-C-C with hydroxyl groups\n'
                 '\n'
                 '2. Improve substituent detection:\n'
                 '- Create more precise SMARTS patterns for acyl/alkyl/alkenyl '
                 'groups\n'
                 '- Better handle the different substitution positions\n'
                 '- Add additional validation checks for chain length\n'
                 '\n'
                 '3. Add molecular property checks:\n'
                 '- Validate molecular weight ranges\n'
                 '- Check carbon count\n'
                 '- Improve oxygen counting logic\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 77,
    'num_false_positives': 27,
    'num_true_negatives': 142171,
    'num_false_negatives': 25,
    'num_negatives': None,
    'precision': 0.7403846153846154,
    'recall': 0.7549019607843137,
    'f1': 0.7475728155339806,
    'accuracy': 0.9996345748418833,
    'negative_predictive_value': 0.9998241863343554}