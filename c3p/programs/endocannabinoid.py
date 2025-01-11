"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:67194 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are characterized by long hydrocarbon chains connected to
    ethanolamide, glycerol, or similar groups via ester, amide, or ether bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible ethanolamide pattern
    ethanolamide_pattern = Chem.MolFromSmarts("[NX3H][CX4H2][CX4H2][OX2H]")
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)

    # Glycerol or similar backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    
    # Ether linkage pattern
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    has_ether = mol.HasSubstructMatch(ether_pattern)

    # Check for at least one of the key functional groups
    if not (has_ethanolamide or has_glycerol or has_ether):
        return False, "No ethanolamide, glycerol, or ether group found"

    # Check for long hydrocarbon chain (at least 12 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for ester, amide, or ether linkage between chain and functional group
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX3](=[OX1])")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern) or has_ether):
        return False, "No ester, amide, or ether linkage between chain and functional group"

    # Check molecular weight (typically >250 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for endocannabinoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12:
        return False, "Too few carbons for endocannabinoid"
    if o_count < 2:
        return False, "Too few oxygens for endocannabinoid"

    # Remove double bond requirement to include saturated endocannabinoids
    # Additional check: at least one functional group oxygen
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No functional group oxygen found"

    return True, "Contains long hydrocarbon chain connected to ethanolamide, glycerol, or similar group via ester, amide, or ether bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:67197',
                          'name': 'endocannabinoid',
                          'definition': 'A class of cannabinoids present in '
                                        'mammalian biological fluids and '
                                        'tissues that activate cannabinoid '
                                        'receptors.',
                          'parents': ['CHEBI:18059', 'CHEBI:67194'],
                          'xrefs': ['Patent:US2012165388'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.17777777777777778 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O=C(O[C@@H]([C@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C '
               'NAME: Roselipin 3D REASON: CORRECT Contains long hydrocarbon '
               'chain connected to ethanolamide, glycerol, or similar group '
               'via ester, amide, or ether bond\n'
               ' * SMILES: C(CCC)C/C=C\\C/C=C\\CC1C(C/C=C\\CCCC(NCCO)=O)O1 '
               'NAME: N-[(5Z,11Z,14Z)-8,9-epoxyicosatrienoyl]ethanolamine '
               'REASON: CORRECT Contains long hydrocarbon chain connected to '
               'ethanolamide, glycerol, or similar group via ester, amide, or '
               'ether bond\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCO NAME: '
               'anandamide REASON: CORRECT Contains long hydrocarbon chain '
               'connected to ethanolamide, glycerol, or similar group via '
               'ester, amide, or ether bond\n'
               ' * SMILES: C(CCCCCCCCCCCCC/C=C\\CCCCCCCC)(NCCO)=O NAME: '
               'N-(15Z-tetracosenoyl)ethanolamine REASON: CORRECT Contains '
               'long hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCOC(CO)CO NAME: '
               '2-arachidonyl glyceryl ether REASON: CORRECT Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@H](O)[C@H](O)CO)[C@H](O)CO)/C(=C/[C@@H]([C@H](O)/C(=C/[C@@H]([C@H](O)/C(=C/[C@@H]([C@H](O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)/C(=C/[C@H](C[C@H](CC)C)C)/C)C)/C)C)/C)C)/C '
               'NAME: TMC-171A REASON: CORRECT Contains long hydrocarbon chain '
               'connected to ethanolamide, glycerol, or similar group via '
               'ester, amide, or ether bond\n'
               ' * SMILES: C(CCCO)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(NCCO)=O '
               'NAME: N-(20-hydroxyarachidonoyl)ethanolamine REASON: CORRECT '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: CCCCCCCC\\C=C/CCCCCCCC(=O)NCCO NAME: oleoyl '
               'ethanolamide REASON: CORRECT Contains long hydrocarbon chain '
               'connected to ethanolamide, glycerol, or similar group via '
               'ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C '
               'NAME: Roselipin 3E REASON: CORRECT Contains long hydrocarbon '
               'chain connected to ethanolamide, glycerol, or similar group '
               'via ester, amide, or ether bond\n'
               ' * SMILES: C(CCC)CC1C(C/C=C\\C/C=C\\C/C=C\\CCCC(NCCO)=O)O1 '
               'NAME: N-[(5Z,8Z,11Z)-14,15-epoxyicosatrienoyl]ethanolamine '
               'REASON: CORRECT Contains long hydrocarbon chain connected to '
               'ethanolamide, glycerol, or similar group via ester, amide, or '
               'ether bond\n'
               ' * SMILES: C(CCC)C/C=C\\C/C=C\\C/C=C\\CC1C(CCCC(NCCO)=O)O1 '
               'NAME: N-[(8Z,11Z,14Z)-5,6-epoxyicosatrienoyl]ethanolamine '
               'REASON: CORRECT Contains long hydrocarbon chain connected to '
               'ethanolamide, glycerol, or similar group via ester, amide, or '
               'ether bond\n'
               ' * SMILES: CCCCCC\\C=C/CCCCCCCC(=O)NCCO NAME: palmitoleoyl '
               'ethanolamide REASON: CORRECT Contains long hydrocarbon chain '
               'connected to ethanolamide, glycerol, or similar group via '
               'ester, amide, or ether bond\n'
               'False positives: SMILES: '
               'OC1C(C2(C(C3C(C(=O)C2)(C4C(=CC3)C(C(=O)C(OC5OC(C(O)C(O)C5O)CO)=C4)(C)C)C)(C1)C)C)C(O)(C)C(=O)/C=C/C(O)(C)C '
               'NAME: Cucurbitacin I 2-glucoside REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1NC=2C(=O)C(=C(OC)C(C2)=O)[C@H](O)[C@H](C[C@H](OC)[C@H](O)[C@H](C=C([C@@H]([C@H](C=CC=C1C)OC)OC(=O)N)C)C)C '
               'NAME: 15-Hydroxygeldanamycin REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               '[H][C@@]12C\\C=C(C)\\[C@@H](O[C@H]3C[C@H](OC)[C@@H](O[C@H]4C[C@H](OC)[C@@H](NC)[C@H](C)O4)[C@H](C)O3)[C@@H](C)\\C=C\\C=C3/CO[C@]4([H])[C@H](O)C(C)=C[C@@]([H])(C(=O)O[C@@H](C1)C[C@]1(O2)O[C@]([H])([C@@H](C)CC)[C@@H](C)C=C1)[C@]34O '
               'NAME: emamectin B1a REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O[C@@H]1[C@]([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@](CC4)(CC[C@H]([C@]5(O)C)C)C(O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)=O)[H])=CC3)C)(CC2)C)[H])(CC1)C)[H])(CO)C '
               'NAME: Pedunculoside REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C2=C1[C@@H](O)[C@@H]3[C@@]4([C@H](O[C@@H](C(O)(C)C)CC4)CC[C@]3(O2)C)C)[C@@H](C(=O)[C@H](CC)C)C '
               'NAME: Cochlioquinone K REASON: WRONGLY CLASSIFIED Contains '
               'long hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'CO[C@@H]1C[C@H](C[C@@H](C)[C@@H]2CC(=O)[C@H](C)\\C=C(C)\\[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\\C=C\\C=C\\C=C(C)\\[C@H](C[C@@H]3CC[C@@H](C)[C@@](O)(O3)C(=O)C(=O)N3CCCC[C@H]3C(=O)O2)OC)CC[C@H]1OC(=O)C(C)(CO)CO '
               'NAME: Temsirolimus REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: O1C(C1CC/C=C/CCCCCCCCCCCCCC)CCCCCCCCCCC2=CC(OC2=O)C '
               'NAME: Epomuricenin B REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (5Z)-tetradecenoyl-CoA REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C(O)C=CC=C[C@H]1C(=C[C@@H](OC)[C@@H]2[C@@H]1[C@@H](C[C@@H](C2)CO)C)C '
               'NAME: Tanzawaic acid S REASON: WRONGLY CLASSIFIED Contains '
               'long hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1OC(C(OC)C=CC=C(CC(C(C(C=C(C=C1OC)C)C)O)C)C)C(C(O)C(C2(OC(C=CC=CC)C(C)C(C2)OC3OC(C(O)C(C3)O)C)O)C)C '
               'NAME: PC-766b REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C(O[C@H]1OC[C@@]23[C@]14[C@]([C@H]5CC(C)(C)C[C@H]5[C@@H]2O3)(C)C4)CCCCCCC/C=C\\CCCCCCCC '
               'NAME: Oleoyl velutinal REASON: WRONGLY CLASSIFIED Contains '
               'long hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1O[C@](OC)(C2=C(O)C=CC(=C2)O)C=C1CC/C=C(/CC/C=C(/C(=O)O)\\C)\\C '
               'NAME: (-)-ganotheaecolumol G REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C(O)/C(=C/CC(OC(=O)C)[C@H]([C@@H]1[C@@]2([C@@](C3=C([C@@]4([C@H](C([C@H](OC(=O)C)CC4)(C)C)C[C@H]3OCC)C)CC2)(C)[C@H](C1)OC(=O)C)C)C)/C '
               'NAME: 7-O-ethyl ganoderic acid O REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain connected to ethanolamide, '
               'glycerol, or similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O1C2C(CC(CC=C(C)C)=C(C)C)(C(O)=C(C(=O)C2C(=O)C(C1C)C)C(O)(C3=CC=CC=C3)C(=O)C)CC=C(C)C '
               'NAME: Mesuafenol REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'CC(\\C=C\\[C@]1(O)[C@@]2(C)CO[C@]1(C)C[C@H](C2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)=C\\C(O)=O '
               'NAME: dihydrophaseic acid 4-O-beta-D-glucoside REASON: WRONGLY '
               'CLASSIFIED Contains long hydrocarbon chain connected to '
               'ethanolamide, glycerol, or similar group via ester, amide, or '
               'ether bond\n'
               ' * SMILES: '
               'P(OC[C@H](O)COCC(OC)CCC/C=C\\CCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: 1-(2-methoxy-6Z-hexadecenyl)-sn-glycero-3-phosphoserine '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain '
               'connected to ethanolamide, glycerol, or similar group via '
               'ester, amide, or ether bond\n'
               ' * SMILES: '
               '[H][C@]1(CSC(=O)N1)[C@@]1(O)C[C@H]2C[C@@H](CC[C@H](C)\\C=C/C=C/CC\\C(C)=C/C(=O)O2)O1 '
               'NAME: latrunculin A REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1N([C@@](C(=O)O)(CO)[C@]([C@H]1C)(O)[C@@H](OC)C[C@H]([C@@H](O)C=CC=CCNC(=O)C([C@H](O)C(=CC=CC=CCC=2OC=NC2)C)(C)C)C)C '
               'NAME: Oxazolomycin A2 REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1C=CC(=CC(C(OC(=O)CC(C(C(C(CC1C)C)OC2OC(CC(C2O)N(C)C)C)C)O)CC)CO)C '
               'NAME: Izenamicin B2 REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               '[C@]12([C@@H]([C@](O[C@]3(O1)[C@H]([C@H](O[C@]3(CC2=O)C)C)C([O-])=O)([C@H](C)/C=C(\\C)/C=C/C(=C/4\\C(NCC4=O)=O)/[O-])[H])C)[H] '
               'NAME: nocamycin E(2-) REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               '[H][C@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)[C@]([H])([C@H](C[C@]34O)OC(C)=O)C3=CC(=O)OC3)[C@@]1(C)CC[C@@H](C2)O[C@H]1C[C@H](O)[C@H](O[C@H]2C[C@H](O)[C@H](O[C@H]3C[C@H](O)[C@H](O)[C@@H](C)O3)[C@@H](C)O2)[C@@H](C)O1 '
               'NAME: 16-O-acetylgitoxin REASON: WRONGLY CLASSIFIED Contains '
               'long hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1O[C@@H](CC=CC=CC(OC2OC(C(N(C)C)CC2)C)[C@@H](C[C@@H]([C@@H]([C@H](C(C1)O)OC)OC3OC(C(OC4OC(C(OC(=O)CC(C)C)C(C4)(O)C)C)C(C3O)N(C)C)C)CC=O)C)C '
               'NAME: Shengjimycin E REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1/C(/OC(=C1)C)=C/C(C(O)C(C(OC)CCC(O)(/C=C(/C2OC(C(C(=O)O)C)CC(=C2)CC)\\C)C)C)C '
               'NAME: Tuscoron A REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O1C(C2C3(OC4=C(C(=O)C3=CCC2)C(O)=C5C(OC(C=C5)(CCC=C(C)C)C)=C4CC=C(C)C)[C@]1(C/C=C(\\C)/C(O)=O)C(O)=O)(C)C '
               'NAME: Garcinolic acid REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               ' * SMILES: '
               'O=C1OC=2C(=O)O[C@]3(C2O)[C@@](C=C(C(=O)O)[C@@H](C3)CC)(CCCC[C@H]4[C@@]1([C@H]5[C@@H]([C@@H](O[C@@H]6OC([C@H](OC(=O)C7=C(OC)C=CC=C7C)[C@@H]([C@@H]6O)O[C@@H]8O[C@@H]([C@H](OC)[C@@H](C8)O)C)C)[C@@H](O)CC5)C=C4C)C)C '
               'NAME: PA-46101 A REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain connected to ethanolamide, glycerol, or '
               'similar group via ester, amide, or ether bond\n'
               'False negatives: SMILES: CCCCCCCCCCCCCCCC(=O)NCCO NAME: '
               'palmitoyl ethanolamide REASON: MISSED No double bonds in '
               'hydrocarbon chain\n'
               ' * SMILES: O(C(CO)CO)C(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC '
               'NAME: MG(0:0/20:3(5Z,8Z,11Z)/0:0) REASON: MISSED No '
               'ethanolamide, glycerol, or ether group found\n'
               ' * SMILES: '
               'O[C@@H](CCCCC)/C=C/C=C\\C/C=C\\C/C=C\\CCCC(OC(CO)CO)=O NAME: '
               '15-HETE-G REASON: MISSED No ethanolamide, glycerol, or ether '
               'group found\n'
               ' * SMILES: '
               '[H]C(CO)(CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC NAME: '
               '2-arachidonoylglycerol REASON: MISSED No ethanolamide, '
               'glycerol, or ether group found\n'
               ' * SMILES: '
               'O(C(C(O)([2H])[2H])(C(O)([2H])[2H])[2H])C(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC '
               'NAME: 2-arachidonoyl glycerol-d5 REASON: MISSED No '
               'ethanolamide, glycerol, or ether group found\n'
               ' * SMILES: CCCCCC=CCC=CCC=CCC=CCCCC(=O)OC(CO)CO NAME: '
               'eicosa-5,8,11,14-tetraenoic acid 1,3-dihydroxypropan-2-yl '
               'ester REASON: MISSED No ethanolamide, glycerol, or ether group '
               'found\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO NAME: '
               'MG(0:0/20:2(11Z,14Z)/0:0) REASON: MISSED No ethanolamide, '
               'glycerol, or ether group found\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\CCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/20:1(11Z)/0:0) REASON: MISSED No ethanolamide, '
               'glycerol, or ether group found\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:3(11Z,14Z,17Z)/0:0) REASON: MISSED No '
               'ethanolamide, glycerol, or ether group found\n'
               ' * SMILES: O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:4(8Z,11Z,14Z,17Z)/0:0) REASON: MISSED No '
               'ethanolamide, glycerol, or ether group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cl.C=1C=C(OC)C(OC)=C(C1CN2CCNCC2)OC.Cl',
                                     'name': 'Trimetazidine hydrochloride',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'},
                                 {   'smiles': '[H]S(=O)(=O)C[C@H](N)C(O)=O',
                                     'name': 'L-cysteine-S-dioxide',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C=C(CO[C@@H]1CN(C)CC3=CC=CC=C3OC)N=N2)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'No long hydrocarbon chain '
                                               'found'},
                                 {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'},
                                 {   'smiles': 'O=C1[C@H]([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@H](O)CC(=C)C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)C',
                                     'name': 'Gilvsin C',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'},
                                 {   'smiles': 'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C(=O)C',
                                     'name': '1-[(4S,7S,8S)-5-acetyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-phenylurea',
                                     'reason': 'No long hydrocarbon chain '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3=NC=CN=C3',
                                     'name': 'N-[(4S,7R,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[oxo(2-pyrazinyl)methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No long hydrocarbon chain '
                                               'found'},
                                 {   'smiles': 'COc1ccc(cc1)\\C=C(NC=O)\\C(NC=O)=C\\c1ccc(OS(O)(=O)=O)cc1',
                                     'name': 'Fumiformamide',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'},
                                 {   'smiles': 'CCCCC(CC)COCCCN',
                                     'name': '3-(2-ethylhexoxy)propan-1-amine',
                                     'reason': 'No long hydrocarbon chain '
                                               'found'},
                                 {   'smiles': 'N(=CC1=CC=CC=C1)CCC2=CC=CC=C2',
                                     'name': 'N-benzylidene-2-phenylethanamine',
                                     'reason': 'No ethanolamide, glycerol, or '
                                               'ether group found'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'The previous program had several issues that led to a low F1 '
                 'score:\n'
                 '\n'
                 '1. **Overly strict double bond requirement**: The program '
                 'required at least one double bond in the hydrocarbon chain, '
                 'which excluded valid endocannabinoids like palmitoyl '
                 'ethanolamide that have fully saturated chains.\n'
                 '\n'
                 '2. **Incomplete functional group detection**: The program '
                 "missed some valid endocannabinoids because it didn't "
                 'properly detect all glycerol-containing structures or ether '
                 'linkages.\n'
                 '\n'
                 '3. **Incorrect pattern matching**: The ethanolamide pattern '
                 'was too specific and missed some valid structures.\n'
                 '\n'
                 '4. **Insufficient structural flexibility**: The program '
                 "didn't account for the full diversity of endocannabinoid "
                 'structures, particularly in how the functional groups are '
                 'connected to the hydrocarbon chain.\n'
                 '\n'
                 "Here's an improved version that addresses these issues:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 23,
    'num_false_positives': 100,
    'num_true_negatives': 619,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.18699186991869918,
    'recall': 1.0,
    'f1': 0.3150684931506849,
    'accuracy': 0.8652291105121294,
    'negative_predictive_value': 1.0}