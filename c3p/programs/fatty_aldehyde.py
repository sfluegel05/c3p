"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde with a carbonyl group at one end of an aliphatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal aldehyde group [C]=O pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"
    
    # Check if surrounds are aliphatic and chainlike
    for match in aldehyde_matches:
        start_atom_idx = match[0]  # Assuming terminal bond to aldehyde is [C]=O
        visited = set()
        carbon_chain = 0

        def is_aliphatic(atom):
            return atom.GetAtomicNum() == 6 and not atom.IsInRing()

        def traverse(atom_idx):
            nonlocal carbon_chain
            if atom_idx not in visited:
                visited.add(atom_idx)
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if is_aliphatic(neighbor):
                        carbon_chain += 1
                        traverse(neighbor.GetIdx())

        traverse(start_atom_idx)

        # Adjust chain length criteria
        if carbon_chain >= 5:
            return True, "Valid fatty aldehyde: Terminal aldehyde group with a suitable aliphatic chain"
    
    return False, f"Carbon chain too short for typical fatty aldehyde (found {carbon_chain} carbons)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35746',
                          'name': 'fatty aldehyde',
                          'definition': 'An aldehyde formally arising from '
                                        'reduction of the carboxylic acid '
                                        'group of its corresponding fatty '
                                        'acid, having a carbonyl group at one '
                                        'end of the carbon chain.',
                          'parents': ['CHEBI:59768'],
                          'xrefs': [   'PMID:15210368',
                                       'PMID:21341652',
                                       'PMID:21347727'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7961165048543689 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: O=C\\C=C\\CCCCCCCCC/C=C\\CCCC NAME: '
               '2E,13Z-Octadecadienal REASON: CORRECT Valid fatty aldehyde: '
               'Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: O=CCCCCCCC/C=C\\C=C/CCCC NAME: '
               '9Z,11Z-Hexadecadienal REASON: CORRECT Valid fatty aldehyde: '
               'Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: OC(CCCC=O)C NAME: 5-hydroxy caproaldehyde REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: CCCC/C=C/C=O NAME: (E)-hept-2-enal REASON: CORRECT '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: O=CCCCCC/C=C\\CCCCCC NAME: 7Z-Tetradecenal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: O=CCC/C=C\\CCCCCCC NAME: (Z)-4-Dodecenal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: C(/C=C/C=C/C)([H])=O NAME: (E,E)-2,4-hexadienal '
               'REASON: CORRECT Valid fatty aldehyde: Terminal aldehyde group '
               'with a long aliphatic chain\n'
               ' * SMILES: O=C\\C=C\\CC/C=C\\C NAME: 2E,6Z-Octadienal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: C(CCCCCCCC)/C=C/C=O NAME: (E)-dodec-2-enal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: O=CCCCCC=O NAME: Adipaldehyde REASON: CORRECT Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: [H]C(=O)CCCCCCCCCCCC NAME: tridecanal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: O=CC(CCCCCCCCC)C NAME: 2-methyl-Undecanal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: O=CCCCCCCCCC/C=C/C=C/CC NAME: 11,13-hexadecadienal '
               'REASON: CORRECT Valid fatty aldehyde: Terminal aldehyde group '
               'with a long aliphatic chain\n'
               ' * SMILES: O=C/C(=C/CCCCCC/C=C\\CCCCCCCC)/CCCCC/C=C\\CCCCCCCC '
               'NAME: (2E,10Z)-2-[(Z)-pentadec-6-enyl]nonadeca-2,10-dienal '
               'REASON: CORRECT Valid fatty aldehyde: Terminal aldehyde group '
               'with a long aliphatic chain\n'
               ' * SMILES: O=CC(CCCCCC)C NAME: 2-Methyloctanal REASON: CORRECT '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: C(CCCCCCC)C(=CC=O)[H] NAME: 2-undecenal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: O=CC/C=C\\CCCCC NAME: (3Z)-nonenal REASON: CORRECT '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: O=C/C=C\\C=C\\C=O NAME: (E,E)-2,4-Hexadienedial '
               'REASON: CORRECT Valid fatty aldehyde: Terminal aldehyde group '
               'with a long aliphatic chain\n'
               ' * SMILES: O=CCCCC/C=C/C NAME: 6-octenal REASON: CORRECT Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: [H]C(=O)C([H])=C([H])C([H])=C([H])C([H])=C([H])CC '
               'NAME: nona-2,4,6-trienal REASON: CORRECT Valid fatty aldehyde: '
               'Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: O=CCC(CC(C)C)C NAME: 3,5-Dimethylhexanal REASON: '
               'CORRECT Valid fatty aldehyde: Terminal aldehyde group with a '
               'long aliphatic chain\n'
               ' * SMILES: '
               'O=C/C(=C/CCCCCC/C=C\\C/C=C\\CCCCC)/CCCCC/C=C\\CCCCCCCC NAME: '
               '(2E,10Z,13Z)-2-[(Z)-pentadec-6-enyl]nonadeca-2,10,13-trienal '
               'REASON: CORRECT Valid fatty aldehyde: Terminal aldehyde group '
               'with a long aliphatic chain\n'
               ' * SMILES: O=CCC/C=C\\CC NAME: (Z)-4-Heptenal REASON: CORRECT '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: O=C/C=C/C=C/C=C/C=C/C=C/C=C/C/C=C\\C NAME: '
               'Thailandene B REASON: CORRECT Valid fatty aldehyde: Terminal '
               'aldehyde group with a long aliphatic chain\n'
               ' * SMILES: O=C/C=C/C=C/C=C/C=C/C=C/C=C/C=O NAME: '
               '2,4,6,8,10,12-tetradeca-1,14-hexenedial REASON: CORRECT Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               'False positives: SMILES: '
               '[NH3+]C[C@H]1O[C@H](O[C@@H]2[C@@H]([NH3+])C[C@@H]([NH3+])[C@H](O)[C@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@H]3O[C@H](C=O)[C@@H](O)[C@H](O)[C@H]3[NH3+])[C@H]2O)[C@H]([NH3+])[C@@H](O)[C@@H]1O '
               "NAME: 6'''-oxoneomycin C(5+) REASON: WRONGLY CLASSIFIED Valid "
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: '
               '[H]C(=O)[C@H](N)[C@@H](O)[C@H](O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1N)[C@H](O)CO '
               'NAME: beta-D-glucosaminyl-(1->4)-aldehydo-D-glucosamine '
               'REASON: WRONGLY CLASSIFIED Valid fatty aldehyde: Terminal '
               'aldehyde group with a long aliphatic chain\n'
               ' * SMILES: N[C@@H](CC[C@@H](O)C=O)C(O)=O NAME: '
               'erythro-L-hydroxyallysine REASON: WRONGLY CLASSIFIED Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: O=C1O/C(=C/C=O)/C(=C1)N NAME: Basidalin REASON: '
               'WRONGLY CLASSIFIED Valid fatty aldehyde: Terminal aldehyde '
               'group with a long aliphatic chain\n'
               ' * SMILES: '
               'C([C@H]([C@H]([C@@H]([C@@H](COP(=O)([O-])[O-])O)O)O)NC(=O)C)=O '
               'NAME: N-acetyl-D-mannosamine 6-phosphate(2-) REASON: WRONGLY '
               'CLASSIFIED Valid fatty aldehyde: Terminal aldehyde group with '
               'a long aliphatic chain\n'
               ' * SMILES: '
               'OC12C3C(C4(C(O)(CC3)CC(OC5OC(C(OC6OC(C(O)C(O)C6O)COC7OC(C(O)C(O)C7O)CO)C(O)C5)C)CC4)C=O)CCC1(C(CC2)C=8COC(=O)C8)C '
               'NAME: Olitoriusin REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'O1C2(C1(CC(O)CC2(C)C)C)\\C=C\\C(\\C)=C/C=C/C(/C)=C/C=C/C=C(\\C)/C=C/C=C(\\C)/C=O '
               'NAME: Beta-Citraurin epoxide REASON: WRONGLY CLASSIFIED Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: '
               'O[C@@H]1[C@]2([C@](C2)(C(=C3[C@]1(CC(C3)(C)C)[H])C=O)C=O)C '
               'NAME: (-)-Merulidial REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'O=C(O[C@@H]1[C@H](OC)[C@H](O[C@H]([C@H]1O)OC[C@]23[C@@]4(C(C(C)C)=C[C@@H]2C[C@]4(C=O)[C@@H]5CC[C@H]([C@H]5C3)C)C(=O)O)CO)/C(=C\\C=C\\C)/C '
               'NAME: '
               '(1R,2S,4R,5R,8R,9S,11S)-9-ormyl-2-[[(2R,3S,4S,5R,6R)-3-hydroxy-6-(hydroxymethyl)-5-methoxy-4-[(2Z,4E)-2-methylhexa-2,4-dienoyl]oxyoxan-2-yl]oxymethyl]-5-methyl-13-propan-2-yltetracyclo[7.4.0.02,11.04,8]tridec-12-ene-1-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Valid fatty aldehyde: Terminal '
               'aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               '[H][C@]1(OC[C@@H](O)[C@H](O)[C@H]1O)O[C@H]1[C@H](O)[C@H](O[C@@]([H])(O[C@H]2CC[C@@]3(C)[C@@]([H])(CC[C@]4(C)[C@]3([H])CC=C3[C@]5([H])CC(C)(C)CC[C@@]5([C@H](O)C[C@@]43C)C(=O)O[C@@H]3O[C@H](C)[C@H](O)[C@H](O)[C@H]3O[C@]3([H])O[C@@H](C)[C@H](O[C@]4([H])OC[C@@H](O)[C@H](O[C@]5([H])OC[C@H](O)[C@H](O)[C@H]5O[C@]5([H])OC[C@H](O)[C@H](O)[C@H]5O)[C@H]4O)[C@@H](O)[C@H]3O)[C@]2(C)C=O)[C@@H]1O[C@]1([H])O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)C(O)=O '
               'NAME: gypsosaponin A REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCC(=O)/C=C/C=O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PKOOA-PG REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCC(=O)/C=C/C=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(O)(O)=O '
               'NAME: OKOOA-PA REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: [H]C(=O)C(\\CC([O-])=O)=C/C=C(\\O)C([O-])=O NAME: '
               '5-formyl-2-hydroxyhepta-2,4-dienedioate REASON: WRONGLY '
               'CLASSIFIED Valid fatty aldehyde: Terminal aldehyde group with '
               'a long aliphatic chain\n'
               ' * SMILES: '
               'CC[C@H]1OC(=O)C[C@@H](O)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)[C@@H](O)[C@@H]([C@H]2O)[NH+](C)C)[C@@H](CC=O)C[C@@H](C)C(=O)\\C=C\\C(C)=C\\[C@@H]1CO '
               'NAME: 5-O-mycaminosyltylonolide(1+) REASON: WRONGLY CLASSIFIED '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: '
               'O=C(O)C12C(=CC3C1(C[C@H]4[C@H](C2(C=O)C3)CC[C@@]4(O)C)CO)C(C)C '
               'NAME: 7-Hydroxysordaricin REASON: WRONGLY CLASSIFIED Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: '
               '[H][C@@](C)(O)[C@]([H])(N)[C@@](C)(O)[C@@]([H])(OC)C=O NAME: '
               'kansosamine REASON: WRONGLY CLASSIFIED Valid fatty aldehyde: '
               'Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'C[C@H](CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C)C=O '
               'NAME: (25R)-3beta-hydroxycholest-5-en-26-al REASON: WRONGLY '
               'CLASSIFIED Valid fatty aldehyde: Terminal aldehyde group with '
               'a long aliphatic chain\n'
               ' * SMILES: '
               'O=C1OC(/C=C/C(=C\\C(C/C=C/C(=C/C(C(=O)C(C2OC(CC=O)(C)CC2C)C)C)/C)C)/CC)C(C)C=C1 '
               'NAME: Leptofuranin D REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'C(/C=C\\CCCCC=O)\\C=C/C=C/C=C/[C@H]([C@@H](O)CCCC(=O)[O-])SC[C@H]([NH3+])C(=O)[O-] '
               'NAME: 20-oxoleukotriene E4(1-) REASON: WRONGLY CLASSIFIED '
               'Valid fatty aldehyde: Terminal aldehyde group with a long '
               'aliphatic chain\n'
               ' * SMILES: '
               'O=C1O[C@H](O)[C@@]23[C@@]1([C@@H]4[C@@](O)(CC(C4)(C)C)C=C2C=O)C3 '
               'NAME: 9-Hydroxymarasmic acid REASON: WRONGLY CLASSIFIED Valid '
               'fatty aldehyde: Terminal aldehyde group with a long aliphatic '
               'chain\n'
               ' * SMILES: '
               'O=C1OC(CC(O)C=CC=CC(CC(C(C(C(C1)OC(=O)CC)OC)OC2O[C@@H]([C@@H](OC3O[C@H]([C@H](OC(=O)CCC)[C@@](C3)(O)C)C)[C@@H]([C@H]2O)N(C)C)C)CC=O)C)C '
               'NAME: SF-837 A2 REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: [H]C(=O)[C@H](O)C[C@H](O)[C@@H](C)O NAME: '
               'aldehydo-D-paratose REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'O=C/C(=C/C=C/[C@H]([C@@H]1[C@@]2([C@@](C3=C([C@@]4([C@H](C([C@H](O)CC4)(C)C)CC3)C)CC2)(C)CC1)C)C)/C '
               'NAME: Astraeusin C REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               '[Fe].O=CC1=CC(=NO)C(O)C=C1.O=CC1=CC(=NO)C(O)C=C1.ON=C1C=C(C=CC1O)CO '
               'NAME: Viridomycin F REASON: WRONGLY CLASSIFIED Valid fatty '
               'aldehyde: Terminal aldehyde group with a long aliphatic chain\n'
               ' * SMILES: '
               'C1=CC(C([C@]2([C@]1([C@@]3([C@@]([C@@H](C2)OC(C)=O)(C=4[C@@](CC3)([C@@](CC4)(C(CO)CC=O)[H])C)C)[H])C)[H])(C)C)=O '
               'NAME: '
               '(1S,3bR,4R,5aR,9aR,9bR,11aS)-1-(1-hydroxy-4-oxobutan-2-yl)-3b,6,6,9a,11a-pentamethyl-7-oxo-1H,2H,3bH,4H,5H,5aH,6H,7H,9aH,9bH,10H,11H,11aH-cyclopenta[a]phenanthren-4-yl '
               'acetate REASON: WRONGLY CLASSIFIED Valid fatty aldehyde: '
               'Terminal aldehyde group with a long aliphatic chain\n'
               'False negatives: SMILES: CCCCC=O NAME: pentanal REASON: MISSED '
               'Carbon chain too short for typical fatty aldehyde (found 5 '
               'carbons)\n'
               ' * SMILES: CCCC(C)C=O NAME: 2-Methylpentanal REASON: MISSED '
               'Carbon chain too short for typical fatty aldehyde (found 5 '
               'carbons)\n'
               ' * SMILES: [H]C(C(C)(C)C)=O NAME: 2,2-dimethylpropanal REASON: '
               'MISSED Carbon chain too short for typical fatty aldehyde '
               '(found 3 carbons)\n'
               ' * SMILES: [H]C(=O)C(C)C NAME: isobutyraldehyde REASON: MISSED '
               'Carbon chain too short for typical fatty aldehyde (found 3 '
               'carbons)\n'
               ' * SMILES: CC[C@H](C)C(=O)[H] NAME: (S)-2-methylbutanal '
               'REASON: MISSED Carbon chain too short for typical fatty '
               'aldehyde (found 4 carbons)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No terminal aldehyde group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No terminal aldehyde group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Upon analyzing the previous code and the outcomes, the main '
                 'issue appears to be the strict requirement for a "long '
                 'aliphatic carbon chain". The threshold set for fatty '
                 'aldehydes at a minimum of 6 carbon atoms might be too high '
                 'and causing high false-negative rates, especially since '
                 "shorter aldehydes like pentanal are being missed. It's "
                 'important to account for the range of fatty aldehydes that '
                 'may have shorter chains. Additionally, an inappropriate '
                 'constraint might be the avoidance of structures with '
                 'aromatic rings, which may not always be relevant for '
                 'classifying fatty aldehydes.\n'
                 '\n'
                 'To improve the classification:\n'
                 '1. Lower the minimum carbon chain length requirement to '
                 'better capture shorter fatty aldehydes.\n'
                 '2. Remove the constraint on avoiding aromatic rings since it '
                 'might not always rule out a fatty aldehyde correctly.\n'
                 '3. Conduct a more thorough examination of the molecular '
                 'structure to ensure it fits the definition of fatty '
                 'aldehydes and avoid misclassification.\n'
                 '\n'
                 "Here's how to adjust the program to address these points:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 210,
    'num_false_positives': 100,
    'num_true_negatives': 19535,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6774193548387096,
    'recall': 1.0,
    'f1': 0.8076923076923077,
    'accuracy': 0.9949609473418998,
    'negative_predictive_value': 1.0}