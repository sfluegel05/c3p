"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has an L-carnitine backbone with an acyl group attached via an ester bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic carnitine backbone structure:
    # [C@H] or [C@@H] central carbon with:
    # - CH2 with carboxylate
    # - CH2 with trimethylammonium
    # - Oxygen (part of ester)
    carnitine_pattern = Chem.MolFromSmarts('[C@H,C@@H]([CH2][C]([O-])=O)([CH2][N+](C)(C)C)[OX2]')
    
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "Missing carnitine backbone structure"
        
    # Check for ester linkage
    # The oxygen from carnitine should be connected to a C(=O)R group
    ester_pattern = Chem.MolFromSmarts('[C@H,C@@H]([CH2][C]([O-])=O)([CH2][N+])O[C](=O)')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing or incorrect ester linkage"

    # Check for trimethylammonium group
    # Allow for deuterated methyl groups ([CH3,CD3])
    trimethyl_pattern = Chem.MolFromSmarts('[N+]([CH3,CD3])([CH3,CD3])([CH3,CD3])')
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Missing trimethylammonium group"
    
    # Verify charges
    pos_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    if pos_charge != 1 or neg_charge != -1:
        return False, f"Incorrect charge distribution: +{pos_charge}, {neg_charge}"

    # Check stereochemistry
    # Get all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "Missing required stereocenter"
    
    # Find the central carbon (the one with OH, CH2COO-, and CH2N+)
    central_matches = mol.GetSubstructMatches(carnitine_pattern)
    if not central_matches:
        return False, "Cannot determine stereochemistry"
        
    # The central carbon should have R configuration for L-carnitine
    # (Note: The R configuration gives L-carnitine due to CIP priority rules)
    central_carbon_idx = central_matches[0][0]
    found_correct_config = False
    
    for idx, config in chiral_centers:
        if idx == central_carbon_idx:
            if config == 'R':  # R configuration corresponds to L-carnitine
                found_correct_config = True
            break
            
    if not found_correct_config:
        return False, "Incorrect stereochemistry - must be L-configuration"

    return True, "Valid O-acyl-L-carnitine structure with correct L-configuration"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75659',
                          'name': 'O-acyl-L-carnitine',
                          'definition': 'An O-acylcarnitine in which the '
                                        'carnitine component has '
                                        'L-configuration.',
                          'parents': ['CHEBI:17387'],
                          'xrefs': [   'KEGG:C02301',
                                       'MetaCyc:O-Acyl-L-Carnitines'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O([C@@H](C[N+](C([2H])([2H])[2H])(C)C)CC([O-])=O)C(=O)CC NAME: '
               'propionyl-L-carnitine-d3 REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: CC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C NAME: '
               'O-acetyl-L-carnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)\\C=C\\C(O)=O '
               'NAME: O-fumaryl-L-carnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: CCC(C)C(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C NAME: '
               '(R)-2-methylbutyrylcarnitine REASON: MISSED Missing '
               'L-carnitine backbone structure\n'
               ' * SMILES: C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)CCCCC(O)=O '
               'NAME: O-adipoyl-L-carnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C '
               'NAME: O-tetradecanoyl-L-carnitine REASON: MISSED Missing '
               'L-carnitine backbone structure\n'
               ' * SMILES: '
               'O(C(=O)CCC/C=C/CCCCCCCC)[C@H](CC([O-])=O)C[N+](C)(C)C NAME: '
               '(5E)-tetradecenoyl-L-carnitine REASON: MISSED Missing '
               'L-carnitine backbone structure\n'
               ' * SMILES: C[N+](C)(C)C[C@H](CC([O-])=O)OC(=O)C=C NAME: '
               'O-propenoyl-D-carnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C NAME: '
               '(R)-oleoylcarnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               ' * SMILES: CCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C NAME: '
               'O-butanoyl-L-carnitine REASON: MISSED Missing L-carnitine '
               'backbone structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC(=O)CCNC(O)=O',
                                     'name': 'N-carboxy-beta-alanine',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)NC3=CC=CC=C3F)C(=O)N(C[C@H]1OC)C)C)CCC(F)(F)F',
                                     'name': '1-(2-fluorophenyl)-3-[(5S,6R,9R)-5-methoxy-3,6,9-trimethyl-2-oxo-8-(3,3,3-trifluoropropyl)-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'O(C1=CC(=C(O)C=C1)C2=C(O)C=CC(=C2)OC)C',
                                     'name': '5,5′-dimethoxybiphenyl-2,2′-diol',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'COC1=CC=C(CCNC[C@H](O)C2=CC=C(O)C=C2)C=C1OC',
                                     'name': 'denopamine',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': '[Zn++].[S-]C(=S)NCCNC([S-])=S',
                                     'name': 'zineb',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCC(C)C)COC(=O)CCCCCCCCC(CC)C)(O)=O)(O)=O',
                                     'name': 'CDP-DG(a-13:0/i-15:0)',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'O=C1OC(O)C2=C1C[C@](O)([C@H]3CC(C[C@H]3[C@@H]2O)(C)C)C',
                                     'name': 'Lactarolide A',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-2,5-Dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-3-yl]acetamide',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'},
                                 {   'smiles': 'COc1ccc(\\C=C/C2CCC=CC2c2ccc(OC)c(OC)c2)cc1OC',
                                     'name': '3-(3,4-dimethoxyphenyl)-4-[(Z)-3,4-dimethoxystyryl]cyclohex-1-ene',
                                     'reason': 'Missing carnitine backbone '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'C[N+](C)(C)C[C@H](CC([O-])=O)OC(=O)C=C',
                                      'name': 'O-propenoyl-D-carnitine',
                                      'reason': 'Incorrect stereochemistry - '
                                                'must be L-configuration'},
                                  {   'smiles': 'O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]',
                                      'name': 'butyryl-L-carnitine-d3',
                                      'reason': 'Missing or incorrect ester '
                                                'linkage'},
                                  {   'smiles': 'C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)CCCCCCCCCCCCCCC([O-])=O',
                                      'name': 'O-hexadecanedioyl-L-carnitine(1-)',
                                      'reason': 'Incorrect charge '
                                                'distribution: +1, -2'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous attempt and propose '
                 'improvements:\n'
                 '\n'
                 '1. The main issue appears to be that the carnitine backbone '
                 'pattern was too restrictive. The pattern '
                 "`'[C@H,C@@H](CC([O-])=O)(CO)C[N+](C)(C)C'` didn't match "
                 'because:\n'
                 '   - It required a specific CO group that might be part of '
                 'the ester\n'
                 '   - The stereochemistry notation might not have matched all '
                 'valid forms\n'
                 '\n'
                 '2. The stereochemistry check was too simplistic - we need to '
                 'specifically verify L-configuration\n'
                 '\n'
                 "3. The ester pattern `'OC(=O)[#6]'` was incorrect - in "
                 'carnitine derivatives, it should be the other way around '
                 '(the oxygen is attached to the carnitine backbone)\n'
                 '\n'
                 "Here's my improved solution:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 35,
    'num_false_positives': 20,
    'num_true_negatives': 142242,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.6363636363636364,
    'recall': 0.9210526315789473,
    'f1': 0.7526881720430109,
    'accuracy': 0.9998383696416022,
    'negative_predictive_value': 0.9999789096277549}