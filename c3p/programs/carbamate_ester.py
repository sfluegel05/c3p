"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is characterized by the presence of the carbamate functional group (-O-C(=O)-N- or an aromatic variant).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhance SMARTS to account for various bonding and aromatic variations
    carbamate_pattern = Chem.MolFromSmarts("[#8][CX3](=[OX1])[NX3,#7]")
    
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group"

    return False, "No carbamate ester functional group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23003',
                          'name': 'carbamate ester',
                          'definition': 'Any ester of carbamic acid or its '
                                        'N-substituted derivatives.',
                          'parents': ['CHEBI:33308'],
                          'xrefs': ['Wikipedia:Carbamate'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.823931623931624 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C(=O)(OCC1=CN=CC=C1)NCC2=CC=C(C=C2)C(=O)NC3=CC=C(C=C3N)F NAME: '
               'pyridin-3-ylmethyl '
               '{4-[(2-amino-4-fluorophenyl)carbamoyl]benzyl}carbamate REASON: '
               'CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: [Sn](CCCC)(CCCC)(CCCC)CN(CCCCl)C(OC(C)(C)C)=O NAME: '
               'tert-Butyl (3-chloropropyl)((tributylstannyl)methyl)carbamate '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: CC(C)(C)OC(=O)N(C)C1=CC=C(C=C1)C2=NC3=CC=CC=C3S2 '
               'NAME: N-[4-(1,3-benzothiazol-2-yl)phenyl]-N-methylcarbamic '
               'acid tert-butyl ester REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@@H](C[C@@H](O)[C@H](Cc1ccccc1)NC(=O)OC(C)(C)C)Cc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C(N)=O '
               'NAME: L-685,458 REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'CO[C@@H]1\\C=C\\C=C(C)\\Cc2cc(OC)c(Cl)c(c2)N(C)C(=O)C[C@H](OC(=O)[C@H](C)N(C)C(=O)CCS)[C@]2(C)O[C@H]2[C@H](C)[C@@H]2C[C@@]1(O)NC(=O)O2 '
               'NAME: mertansine REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'CC(C)[C@H](NC(=O)OCc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C=O NAME: '
               'Z-Val-Phe-H REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: CCOc1ccc(NC(=O)OC(C)C)cc1OCC NAME: diethofencarb '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: '
               'COC(=O)N[C@@H](C(O)Cc1c(C)nc2n(C)c3n(cnc3c(=O)n12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O)C(=O)OC '
               "NAME: beta-hydroxywybutosine 5'-monophosphate REASON: CORRECT "
               'Contains carbamate ester functional group\n'
               ' * SMILES: CCCC(C)(COC(N)=O)COC(=O)NC(C)C NAME: carisoprodol '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: CNC(=O)ON=C(C)C(C)S(C)(=O)=O NAME: butoxycarboxim '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: CCOC(=O)NCCOc1ccc(Oc2ccccc2)cc1 NAME: fenoxycarb '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: COC(=O)NC(=S)Nc1ccccc1NC(=S)NC(=O)OC NAME: '
               'thiophanate-methyl REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'N(C(=O)O[C@H](C)C1=CC=CC=C1Cl)C2=C(C=3C=CC(CSCCC(=O)O)=CC3)ON=C2C '
               'NAME: (R)-KI16425 REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: CCCCC1=CC2=C(NC(NC(=O)OC)=N2)C=C1 NAME: '
               'parbendazole REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'CO[C@@H]1[C@@H](O)[C@H](NC(=O)CC(C)(C)O)[C@@H](C)O[C@@H]1OCCNC(=O)OCc1ccccc1 '
               'NAME: '
               'N-benzyloxycarbonylaminoethyl-4,6-dideoxy-4-(3-hydroxy-3-methylbutanamido)-2-O-methyl-beta-D-glucopyranoside '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)OC(C)(C)C)\\C=C\\[C@H](CC1=CC=CC=C1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCNC(=O)OCC1=CC=CC=C1)C(=O)NC1CC(C)(C)N([O])C(C)(C)C1 '
               'NAME: XJB-5-131 REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: Cl.CN(C)C(=O)Oc1cc(OC(=O)N(C)C)cc(c1)C(O)CNC(C)(C)C '
               'NAME: bambuterol hydrochloride REASON: CORRECT Contains '
               'carbamate ester functional group\n'
               ' * SMILES: CCNC(=O)[C@H](C)OC(=O)Nc1ccccc1 NAME: carbetamide '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: CSC(C)=NOC(=O)NSNC(=O)ON=C(C)SC NAME: thiodicarb '
               'REASON: CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: '
               'NC[C@H]1O[C@H](O[C@@H]2[C@@H](N)C[C@@H](N)[C@H](O[C@H]3O[C@H](COC(N)=O)[C@@H](O)[C@H](N)[C@H]3O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               "NAME: 6''-O-carbamoylkanamycin A REASON: CORRECT Contains "
               'carbamate ester functional group\n'
               ' * SMILES: '
               'CC(=O)Nc1ccc(cc1)N1C(=O)C2[C@H]3CC[C@@](CC3)(NC(=O)OCC(O)=O)C2C1=O '
               'NAME: '
               '{[(1r,7r)-4-(4-acetamidophenyl)-3,5-dioxo-4-azatricyclo[5.2.2.0(2,6)]undec-1-yl]carbamoyloxy}acetic '
               'acid REASON: CORRECT Contains carbamate ester functional '
               'group\n'
               ' * SMILES: '
               'CO[C@@H]1\\C=C\\C=C(C)\\Cc2cc(OC)c(Cl)c(c2)N(C)C(=O)C[C@H](OC(=O)[C@H](C)N(C)C(C)=O)[C@]2(C)O[C@H]2[C@H](C)[C@@H]2C[C@@]1(O)NC(=O)O2 '
               'NAME: maytansine REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               ' * SMILES: C1=CC=CC(=C1)NC(OC(C)C)=O NAME: propham REASON: '
               'CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: O(C1(CCCCC1)CC#C)C(N)=O NAME: hexapropymate REASON: '
               'CORRECT Contains carbamate ester functional group\n'
               ' * SMILES: '
               '[H][C@]12OCC[C@@]1([H])[C@H](CO2)OC(=O)N[C@@H](Cc1ccccc1)[C@H](O)CN(CC(C)C)S(=O)(=O)c1ccc(N)cc1 '
               'NAME: darunavir REASON: CORRECT Contains carbamate ester '
               'functional group\n'
               'False positives: SMILES: '
               'COC1=CC(=CC(=C1OC(=O)N2CCOCC2)OC)C3C(=C(OC4=C3C(=O)CCC4)N)C#N '
               'NAME: 4-morpholinecarboxylic acid '
               '[4-(2-amino-3-cyano-5-oxo-4,6,7,8-tetrahydro-1-benzopyran-4-yl)-2,6-dimethoxyphenyl] '
               'ester REASON: WRONGLY CLASSIFIED Contains carbamate ester '
               'functional group\n'
               ' * SMILES: P(=O)(N1CC1)(N2CC2)NC(OCC3=CC=CC=C3)=O NAME: '
               'benzodepa REASON: WRONGLY CLASSIFIED Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'CC(C)C[C@@H](C(=O)OC(C)(C)C)NC(=O)C1=C(NC=N1)C(=O)NCCCCCNC(=O)OC(C)(C)C '
               'NAME: '
               '(2S)-4-methyl-2-[[[5-[[5-[[(2-methylpropan-2-yl)oxy-oxomethyl]amino]pentylamino]-oxomethyl]-1H-imidazol-4-yl]-oxomethyl]amino]pentanoic '
               'acid tert-butyl ester REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: '
               'ClC(Cl)CCCC1C2=C(O)C=C([C@H](OC(=O)N)[C@H](CCCCC(CCCC)C3=C(C=C([C@@H]([C@H](CCCC1)C)OC(=O)N)C=C3O)O)C)C=C2O '
               'NAME: Carbamidocyclophane C REASON: WRONGLY CLASSIFIED '
               'Contains carbamate ester functional group\n'
               ' * SMILES: [O-]C(=O)CNC(=O)OCc1ccccc1 NAME: '
               'N-benzyloxycarbonylglycinate REASON: WRONGLY CLASSIFIED '
               'Contains carbamate ester functional group\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)OC4=CC=C(C=C4)C)[C@H](C)CO '
               'NAME: LSM-30826 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'S1C(=NC(CN(C(=O)N[C@@H](C(C)C)C(=O)N[C@H](C[C@H](O)[C@@H](NC(OCC=2SC=NC2)=O)CC3=CC=CC=C3)CC4=CC=CC=C4)C)=C1)C(O)(C)C '
               'NAME: Hydroxy Ritonavir REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)OC4=CC=CC=C4)[C@@H](C)CO '
               'NAME: LSM-30755 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'CCc1c(C)c2[nH]c1cc1nc(cc3[nH]c(cc4nc(cc5[nH]c2c(C)c5CC)c(C)c4CCCOC(=O)N(CCOCCOCCOC)CCOCCOCCOC)c(CC)c3CC)c(CCCOC(=O)N(CCOCCOCCOC)CCOCCOCCOC)c1C '
               'NAME: sapphyrin PCI-2050 REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: O1C[C@@H](NC1=O)CC2=CC=CC=C2 NAME: '
               '(S)-4-BENZYL-2-OXAZOLIDINONE REASON: WRONGLY CLASSIFIED '
               'Contains carbamate ester functional group\n'
               ' * SMILES: '
               'COC([C@@H]1C[C@@]23CCC[N+]4(CCC5(C6=CC=CC=C6N(C(=O)OC)[C@]15CC2)[C@]34[H])[O-])=O '
               'NAME: Pleiocarpoline REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: COCC1(C)OC(=O)N(C1=O)c1cc(Cl)cc(Cl)c1 NAME: '
               'Myclozolin REASON: WRONGLY CLASSIFIED Contains carbamate ester '
               'functional group\n'
               ' * SMILES: CNC(=O)Oc1ccc2N(C)[C@H]3ON(C)CC[C@@]3(C)c2c1 NAME: '
               'Eseridine REASON: WRONGLY CLASSIFIED Contains carbamate ester '
               'functional group\n'
               ' * SMILES: '
               'C1CN(CCC1C(C2=CC3=C(C=C2)OCO3)(C4=CC5=C(C=C4)OCO5)O)C(=O)OC6=CC=C(C=C6)[N+](=O)[O-] '
               'NAME: '
               '4-[bis(1,3-benzodioxol-5-yl)-hydroxymethyl]-1-piperidinecarboxylic '
               'acid (4-nitrophenyl) ester REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: CNC(=O)Oc1cc(C)c(N(C)C)c(C)c1 NAME: Mexacarbate '
               'REASON: WRONGLY CLASSIFIED Contains carbamate ester functional '
               'group\n'
               ' * SMILES: CN1C(=O)OC(C)(C)C1=O NAME: trimethadione REASON: '
               'WRONGLY CLASSIFIED Contains carbamate ester functional group\n'
               ' * SMILES: CCOC(=O)NC1=NN=C(S1)C2=CC(=C(C=C2)OC)OC NAME: '
               'N-[5-(3,4-dimethoxyphenyl)-1,3,4-thiadiazol-2-yl]carbamic acid '
               'ethyl ester REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'N1(C(C(C2=CC=C(OC3=CC=CC=C3)C=C2)(C)OC1=O)=O)NC4=CC=CC=C4 '
               'NAME: '
               '5-methyl-5-(4-phenoxyphenyl)-3-(phenylamino)-1,3-oxazolidine-2,4-dione '
               'REASON: WRONGLY CLASSIFIED Contains carbamate ester functional '
               'group\n'
               ' * SMILES: O(C(=O)N[C@@H](CCCCNC(=O)C=C)C(O)=O)CC1=CC=CC=C1 '
               'NAME: N6-Acryloyl-N2-[(benzyloxy)carbonyl]lysine REASON: '
               'WRONGLY CLASSIFIED Contains carbamate ester functional group\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@H]1CN(C)C(=O)OCC(C)C)C4=CC=CC=C4N2C)[C@H](C)CO '
               'NAME: LSM-8684 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'S(=O)(=O)(O[C@@H]1[C@@H](O)[C@]23NC(N)=N[C@H]2[C@H](COC(=O)N)N=C(N3C1)N)O '
               'NAME: 12beta-Deoxydecarbamoyloxygonyautoxin-3 REASON: WRONGLY '
               'CLASSIFIED Contains carbamate ester functional group\n'
               ' * SMILES: '
               'O=C1N2C3=C(O[C@@H](OC)C(C[C@H]2C)=O)C(=C(OC)C(=C3)O)C[C@H](C[C@H](OC)[C@H](O)[C@H](C=C([C@@H]([C@H](C=CC=C1C)OC)OC(=O)N)C)C)C '
               'NAME: Natalamycin A REASON: WRONGLY CLASSIFIED Contains '
               'carbamate ester functional group\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)OC4=CC=C(C=C4)C)[C@@H](C)CO '
               'NAME: LSM-30631 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@H]1CN(C)C(=O)OC(C)C)[C@H](C)CO '
               'NAME: LSM-9858 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               ' * SMILES: '
               'CC(C)(C)OC(=O)N[C@H]1CCCCCC=C[C@@H]2C[C@]2(NC(=O)[C@@H]3C[C@H](CN3C1=O)OC(=O)N4CC5=C(C4)C(=CC=C5)F)CNS(=O)(=O)C6CC6 '
               'NAME: LSM-5352 REASON: WRONGLY CLASSIFIED Contains carbamate '
               'ester functional group\n'
               'False negatives: SMILES: CCOC(=O)n1ccn(C)c1=S NAME: '
               'carbimazole REASON: MISSED No carbamate ester functional group '
               'found\n'
               ' * SMILES: OCn1c2ccc(Cl)cc2oc1=O NAME: '
               '6-chloro-3-(hydroxymethyl)benzoxazolin-2-one REASON: MISSED No '
               'carbamate ester functional group found\n'
               ' * SMILES: CCOP(=S)(OCC)SCn1c2ccc(Cl)cc2oc1=O NAME: phosalone '
               'REASON: MISSED No carbamate ester functional group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No carbamate ester functional '
                                               'group found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No carbamate ester functional '
                                               'group found'}],
    'sample_false_negatives': [   {   'smiles': 'OCn1c2ccc(Cl)cc2oc1=O',
                                      'name': '6-chloro-3-(hydroxymethyl)benzoxazolin-2-one',
                                      'reason': 'No carbamate ester functional '
                                                'group found'},
                                  {   'smiles': 'CCOP(=S)(OCC)SCn1c2ccc(Cl)cc2oc1=O',
                                      'name': 'phosalone',
                                      'reason': 'No carbamate ester functional '
                                                'group found'}],
    'attempt': 2,
    'reasoning': 'In the previous code, we attempted to identify carbamate '
                 'esters by searching for the specific substructure '
                 '`[-O-C(=O)-N-]` within the given SMILES string. This '
                 'approach relies on the assumption that all carbamate esters '
                 'will contain this exact pattern, but as demonstrated by the '
                 'false positives and negatives, this assumption might not '
                 'always hold.\n'
                 '\n'
                 'Issues to address:\n'
                 '1. **Pattern Flexibility**: The SMARTS pattern '
                 '`[O][C](=[O])[N]` may not catch all variations of the '
                 'carbamate ester group. For instance, the fact that some '
                 'atoms in the chain may bear charges or be attached to '
                 'aromatic systems.\n'
                 '   \n'
                 '2. **False Positives**: Some molecules may coincidentally '
                 'have similar structures but do not form part of a carbamate '
                 'ester.\n'
                 '\n'
                 '3. **False Negatives**: Structures that should have been '
                 'identified as carbamate esters were missed. This suggests '
                 'that the SMARTS pattern might need broadening or that '
                 'further spatial arrangements need to be considered.\n'
                 '\n'
                 'To improve the classification, we can:\n'
                 '- Enhance the SMARTS pattern by taking into account common '
                 'variations in organic compounds, likely adding charges or '
                 'aromaticity specifications.\n'
                 '- Add additional constraints or checks for expected chemical '
                 'properties associated with carbamate esters.\n'
                 '  \n'
                 'Below is an improved implementation for the '
                 '`is_carbamate_ester` function:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 242,
    'num_false_positives': 100,
    'num_true_negatives': 14236,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.7076023391812866,
    'recall': 0.9918032786885246,
    'f1': 0.8259385665529011,
    'accuracy': 0.9930041152263375,
    'negative_predictive_value': 0.9998595308329822}