"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal allylic hydroxy group 
of a polyprenol with 1 mol eq. of phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    diphosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_diphosphate = mol.HasSubstructMatch(diphosphate_pattern)
    
    if not (has_phosphate or has_diphosphate):
        return False, "No phosphate/diphosphate group found"

    # Look for isoprene units with flexible matching for cis/trans configurations
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-[C]=[C]-[CH2]-[CH2]")
    isoprene_pattern2 = Chem.MolFromSmarts("[CH3]-[C](-[CH3])=[C]-[CH2]")
    
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    isoprene_matches2 = len(mol.GetSubstructMatches(isoprene_pattern2))
    total_isoprene = max(isoprene_matches, isoprene_matches2)
    
    if total_isoprene < 1:
        return False, "No isoprene units found"
    
    # Check for allylic oxygen-phosphorus connection with flexible pattern
    allylic_phosphate = Chem.MolFromSmarts("[CH2,CH3]-[C]=[C]-[CH2]-[OX2]P")
    if not mol.HasSubstructMatch(allylic_phosphate):
        return False, "No phosphate group connected to allylic position"
        
    # Count carbons to verify chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:  # Minimum 1 isoprene unit (C5)
        return False, "Carbon chain too short for polyprenol"
        
    # Look for characteristic branching pattern of polyprenols
    methyl_branches = Chem.MolFromSmarts("[CH3]-[C]=[C,CH1]")
    branch_count = len(mol.GetSubstructMatches(methyl_branches))
    
    if branch_count < 1:
        return False, "No methyl branches found"

    # Check for phosphate/diphosphate at terminal position
    terminal_phosphate = Chem.MolFromSmarts("[CH2]-[OX2]P(=O)([OX2])[OX2]")
    terminal_diphosphate = Chem.MolFromSmarts("[CH2]-[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    
    if not (mol.HasSubstructMatch(terminal_phosphate) or mol.HasSubstructMatch(terminal_diphosphate)):
        return False, "Phosphate group not at terminal position"

    # Success case - provide detailed reason
    phosphate_type = "diphosphate" if has_diphosphate else "phosphate"
    return True, f"Contains polyprenol chain with {phosphate_type} group at terminal allylic position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16460',
                          'name': 'polyprenol phosphate',
                          'definition': 'A prenol phosphate resulting from the '
                                        'formal condensation of the terminal '
                                        'allylic hydroxy group of a polyprenol '
                                        'with 1 mol eq. of phosphoric acid.',
                          'parents': ['CHEBI:26250', 'CHEBI:26875'],
                          'xrefs': ['KEGG:C01048'],
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
               'C[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](C)O[C@@H]2O[C@@H]2[C@@H](CO)OC(OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)[C@H](NC(C)=O)[C@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O '
               'NAME: '
               'alpha-D-rhamnosyl-(1->2)-alpha-D-rhamnosyl-(1->4)-N-acetyl-D-glucosaminyl '
               'undecaprenyl diphosphate REASON: MISSED Not enough isoprene '
               'units for a polyprenol\n'
               ' * SMILES: '
               'CC(=O)N[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H]1[C@@H](CO)OC(OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)[C@H](NC(C)=O)[C@H]1O '
               'NAME: '
               'N-acetyl-beta-D-mannosaminyl-(1->4)-N-acetyl-D-glucosaminyl '
               'undecaprenyl diphosphate REASON: MISSED Not enough isoprene '
               'units for a polyprenol\n'
               ' * SMILES: '
               'O(C/C=C(\\CC/C=C(/C)\\CC/C=C(\\CCC=C(C)C)/C)/C)P(OP(=O)(O)O)(=O)O '
               'NAME: nerylneryl diphosphate REASON: MISSED Not enough '
               'isoprene units for a polyprenol\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O '
               'NAME: all-trans-heptaprenyl diphosphate REASON: MISSED Not '
               'enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O '
               'NAME: all-trans-pentaprenyl diphosphate REASON: MISSED Not '
               'enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O '
               'NAME: all-trans-hexaprenyl diphosphate REASON: MISSED Not '
               'enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O '
               'NAME: trans,polycis-decaprenyl phosphate REASON: MISSED Not '
               'enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(=O)OP(O)(O)=O '
               'NAME: ditrans,polycis-dodecaprenyl diphosphate REASON: MISSED '
               'Not enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'C[C@@H](NC(=O)[C@@H](C)NC(=O)[C@H](CCCCNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN)NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@H](O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)C(N)=O)C(O)=O '
               'NAME: '
               'undecaprenyldiphospho-N-acetyl-(N-acetylglucosaminyl)muramoyl-L-alanyl-D-isoglutaminyl-(glycyl)5-L-lysyl-D-alanyl-D-alanine '
               'REASON: MISSED Not enough isoprene units for a polyprenol\n'
               ' * SMILES: '
               'C[C@H]1OC(OP(O)(=O)OP(O)(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)[C@H](NC(C)=O)[C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](NC(C)=O)[C@@H]2NC(C)=O)C(O)=O)[C@H]1O '
               'NAME: '
               '2,3-diacetamido-2,3-dideoxy-beta-D-mannosyl-(1->3)-N-acetyl-D-fucosaminyl '
               'undecaprenyl diphosphate REASON: MISSED Not enough isoprene '
               'units for a polyprenol\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No phosphate/diphosphate group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)=CCC\\C(C)=C\\CCC(C)(OP(O)(=O)OP(O)(O)=O)C=C',
                                      'name': 'nerolidyl diphosphate',
                                      'reason': 'No phosphate group connected '
                                                'to allylic position'},
                                  {   'smiles': 'P(O[C@@H]1OC[C@@H](O)[C@H](O)C1O)(OCCC(CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)C)(O)=O',
                                      'name': 'dolichyl-4-D-xylosyl phosphate',
                                      'reason': 'No phosphate group connected '
                                                'to allylic position'},
                                  {   'smiles': 'CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC\\C(C)=C\\COP(O)(O)=O',
                                      'name': 'phytyl phosphate',
                                      'reason': 'No isoprene units found'},
                                  {   'smiles': 'CC(C)CCCC(C)CCCC(C)CCC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O',
                                      'name': '(E)-3,7,11,15-tetramethylhexadec-2-en-1-yl '
                                              'diphosphate',
                                      'reason': 'No isoprene units found'},
                                  {   'smiles': 'CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC\\C(C)=C\\COP(O)(=O)OP(O)(O)=O',
                                      'name': 'phytyl diphosphate',
                                      'reason': 'No isoprene units found'},
                                  {   'smiles': 'CC(CCOP(O)(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C',
                                      'name': 'beta-D-mannosyl '
                                              'C35-phosphodolichol',
                                      'reason': 'No phosphate group connected '
                                                'to allylic position'},
                                  {   'smiles': 'P(O[C@@H]1OC([C@@H](O)[C@H](O)C1O)CO)(OCCC[C@H](CCC[C@H](CCC[C@H](CCC[C@H](CCC[C@H](CCCCCCCCC)C)C)C)C)C)(O)=O',
                                      'name': 'Mannosyl-1beta-phosphomycoketide '
                                              'C34',
                                      'reason': 'No isoprene units found'}],
    'attempt': 1,
    'reasoning': 'Looking at the failed outcomes, I notice that the program '
                 'incorrectly rejected many valid polyprenol phosphates. The '
                 'main issues appear to be:\n'
                 '\n'
                 '1. The isoprene pattern detection was too strict - using '
                 '"C=C-C(C)-C" misses many valid isoprene units because '
                 'polyprenols can have both cis and trans configurations, and '
                 'the SMILES notation uses different symbols for these (\\, /, '
                 'or none).\n'
                 '\n'
                 '2. The minimum requirements were too restrictive - the code '
                 'required at least 2 isoprene units, but some valid compounds '
                 'like geranyl phosphate (2 units) were missed.\n'
                 '\n'
                 '3. The conjugated double bonds check was incorrect - some '
                 'valid polyprenol phosphates can appear to have conjugated '
                 'bonds in their SMILES representation due to the way '
                 'cis/trans bonds are written.\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 104,
    'num_false_positives': 15,
    'num_true_negatives': 142174,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.8739495798319328,
    'recall': 0.9369369369369369,
    'f1': 0.9043478260869565,
    'accuracy': 0.9998453970484891,
    'negative_predictive_value': 0.9999507669801169}