"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid typically has a long carbon chain with multiple cis double bonds,
    a terminal carboxylic acid group, and may include other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing terminal carboxylic acid group (C(=O)O)"
    
    # Check for multiple cis double bonds (excluding aromatic)
    cis_double_bond_pattern = Chem.MolFromSmarts("[C&!a]/C=C/[C&!a]")
    cis_double_bond_matches = mol.GetSubstructMatches(cis_double_bond_pattern)
    if len(cis_double_bond_matches) < 2:
        return False, f"Insufficient number of cis double bonds, found {len(cis_double_bond_matches)}"
    
    # Include check for specific functional groups present in given examples (e.g., phosphocholine)
    phosphocholine_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if mol.HasSubstructMatch(phosphocholine_pattern):
        return True, "Matches essential fatty acid criteria including phosphocholine functionality"

    # Final check for carbon chain length -> broaden range due to variable sizes in examples or remove
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Lower bound from before removed upper bound for broader match
        return False, f"Carbon chain length {c_count} is too short for essential fatty acids"

    return True, "Matches essential fatty acid criteria with multiple cis double bonds and a carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59549',
                          'name': 'essential fatty acid',
                          'definition': 'Any member of the sub-set of '
                                        'polyunsaturated fatty acid for which '
                                        'there is an absolute dietary '
                                        'requirement.',
                          'parents': ['CHEBI:26208'],
                          'xrefs': [   'PMID:1745654',
                                       'PMID:19034351',
                                       'PMID:7609665'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.4044943820224719 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCCCC(O)=O '
               'NAME: (21Z,24Z,27Z,30Z,33Z)-hexatriacontapentaenoic acid '
               'REASON: CORRECT Matches essential fatty acid criteria with '
               'multiple cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCC(O)=O NAME: '
               '(17Z,20Z,23Z,26Z)-dotriacontatetraenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O '
               'NAME: (6Z,9Z,12Z,15Z,18Z)-tetracosapentaenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC(O)=O NAME: '
               'all-cis-8,11,14,17-icosatetraenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'O(C(=O)\\C(=C(\\C(=C(\\C(=C(\\C(=C(\\C(C(C(C(C(C(C(C(C(C(C([H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])([H])[H])/[H])\\[H])/[H])\\[H])/[H])\\[H])/[H])\\[H])[H] '
               'NAME: eicosatetraenoic acid REASON: CORRECT Matches essential '
               'fatty acid criteria with multiple cis double bonds and a '
               'carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCC(O)=O NAME: '
               '(15Z,18Z,21Z,24Z,27Z)-triacontapentaenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCC(O)=O '
               'NAME: (18Z,21Z,24Z,27Z,30Z,33Z)-hexatriacontahexaenoic acid '
               'REASON: CORRECT Matches essential fatty acid criteria with '
               'multiple cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCC(O)=O '
               'NAME: (16Z,19Z,22Z,25Z,28Z,31Z)-tetratriacontahexaenoic acid '
               'REASON: CORRECT Matches essential fatty acid criteria with '
               'multiple cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC(O)=O NAME: '
               '(8Z,11Z,14Z,17Z,20Z,23Z)-hexacosahexaenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCC(O)=O '
               'NAME: (15Z,18Z,21Z,24Z)-triacontatetraenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: C(=C\\C/C=C\\C/C=C\\CCCCCC(O)=O)\\C/C=C\\CCCCC '
               'NAME: all-cis-docosa-7,10,13,16-tetraenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCC(O)=O '
               'NAME: (13Z,16Z,19Z,22Z)-octacosatetraenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCC(O)=O '
               'NAME: (10Z,13Z,16Z,19Z,22Z,25Z)-octacosahexaenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '(19Z,22Z,25Z,28Z)-tetratriacontatetraenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C=C/C=C/CCCCCCC(O)=O NAME: '
               '(8E,10E,12Z)-octadecatrienoic acid REASON: CORRECT Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O NAME: '
               'gamma-linolenic acid REASON: CORRECT Matches essential fatty '
               'acid criteria with multiple cis double bonds and a carboxylic '
               'acid group\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O NAME: '
               'all-cis-octadeca-6,9,12,15-tetraenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCC(O)=O NAME: '
               '(10Z,13Z,16Z,19Z,22Z)-octacosapentaenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: C(CCC(O)=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CC NAME: '
               'all-cis-icosa-11,14,17-trienoic acid REASON: CORRECT Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCC(O)=O NAME: '
               '(9Z,12Z,15Z,18Z)-tetracosatetraenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCC(O)=O NAME: '
               '(13Z,16Z,19Z,22Z,25Z)-octacosapentaenoic acid REASON: CORRECT '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCC(O)=O '
               'NAME: (12Z,15Z,18Z,21Z,24Z,27Z)-triacontahexaenoic acid '
               'REASON: CORRECT Matches essential fatty acid criteria with '
               'multiple cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCC(O)=O '
               'NAME: (11Z,14Z,17Z,20Z)-hexacosatetraenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCC(O)=O '
               'NAME: (17Z,20Z,23Z,26Z,29Z)-dotriacontapentaenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(O)=O '
               'NAME: (7Z,10Z,13Z,16Z,19Z)-docosapentaenoic acid REASON: '
               'CORRECT Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               'False positives: SMILES: '
               'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H] '
               'NAME: Urceolide REASON: WRONGLY CLASSIFIED Matches essential '
               'fatty acid criteria with multiple cis double bonds and a '
               'carboxylic acid group\n'
               ' * SMILES: '
               'CC(C)[C@@]1(O)O[C@H]2C[C@@]3(C)OC(=CC3=O)\\C(C)=C/[C@H]3OC(=O)[C@@]1(C)[C@H]23 '
               'NAME: Eremantholide A REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'C[C@@H]1C[C@H](O)[C@@H]2[C@@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H](C)C=CC2=C1 '
               'NAME: monacolin J REASON: WRONGLY CLASSIFIED Matches essential '
               'fatty acid criteria with multiple cis double bonds and a '
               'carboxylic acid group\n'
               ' * SMILES: '
               'O=C1O[C@@H](CC=2C1=C(O)C3=C(O)C(C4=C(O)C=5C(=O)C=C(C)C(C5C=C4OC)=O)=C(OC)C=C3C2)C '
               'NAME: Xanthoradone B REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'COC(=O)CCc1c(C)c2=CC3=[N+]4C(=Cc5c(C)c(C=C)c6C=C7[N+]8=C(C=c1n2[Fe--]48n56)[C@]1(CCC(=O)O1)[C@@]7(C)O)C(C=C)=C3C '
               'NAME: cis-heme d hydroxychlorin gamma-spirolactone methyl '
               'ester REASON: WRONGLY CLASSIFIED Matches essential fatty acid '
               'criteria with multiple cis double bonds and a carboxylic acid '
               'group\n'
               ' * SMILES: '
               'O=C1O[C@@H](/C=C/CNC(=O)/C=C\\C=N/OC)[C@@H](O)C=CC=CC=CC=2C1=C(O)C=CC2 '
               'NAME: Oximidine II REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C1C2=C(C[C@]3(OC(=O)C=C3)OC2)C[C@H]4[C@]1(OC(=O)[C@@H]4C(O)CCCCCCC)C '
               'NAME: Sequoiamonascin A REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'C12=C(C=CC=C1)NC=3[C@]([C@@]4(/C(/CN(CC4)CCC32)=C\\C)[H])(C(=O)OC)COC(C)=O '
               'NAME: O-acetyl-15alpha-stemmadenine REASON: WRONGLY CLASSIFIED '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O NAME: '
               '1-(11Z,14Z-eicosadienoyl)-glycero-3-phosphate REASON: WRONGLY '
               'CLASSIFIED Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C(O)[C@]1(O)[C@@H]2C=C3CC[C@@H]4[C@@]([C@]3(C2)CC1)(CC[C@H]([C@]4(CO)C)O)C '
               'NAME: Aphidicolin A25 REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: CC1=C(C(=CC=C1)C)OC(=O)CSC2=NNC(=C3C=CC=CC3=O)O2 '
               'NAME: '
               '2-[[2-(6-oxo-1-cyclohexa-2,4-dienylidene)-3H-1,3,4-oxadiazol-5-yl]thio]acetic '
               'acid (2,6-dimethylphenyl) ester REASON: WRONGLY CLASSIFIED '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C1C([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@H](C(=O)O)CCC=C(C)C)C(C4)=O)(C)CC3)C)C[C@@H]2O)(C)CC1)(C)C '
               'NAME: Formipinic acid E REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'O[C@@H]1[C@]2(C(C3=C([C@@]4(C(C([C@H](OC(=O)CC(O)=O)CC4)(C)C)CC3)C)C1)(CC[C@@]2([C@@H](CCC([C@@H](C)C(O)=O)=C)C)[H])C)C '
               'NAME: '
               '(2R,6R)-6-[(3R,10S,12S,13R,17R)-3-(2-Carboxyacetyl)oxy-12-hydroxy-4,4,10,13,14-pentamethyl-2,3,5,6,7,11,12,15,16,17-decahydro-1H-cyclopenta[a]phenanthren-17-yl]-2-methyl-3-methylideneheptanoic '
               'acid REASON: WRONGLY CLASSIFIED Matches essential fatty acid '
               'criteria with multiple cis double bonds and a carboxylic acid '
               'group\n'
               ' * SMILES: S1C(C1C/C=C/C/C=C/C/C=C/CCCC(O)=O)CCCCC NAME: '
               '(5E,8E,11E)-13-[(2S,3S)-3-Pentylthiiran-2-yl]trideca-5,8,11-trienoic '
               'acid REASON: WRONGLY CLASSIFIED Matches essential fatty acid '
               'criteria with multiple cis double bonds and a carboxylic acid '
               'group\n'
               ' * SMILES: '
               'C1=2C=3C=4[C@](N(C[C@@H](C4)C(N[C@@H](CC)CO)=O)C)(CC1=CN(C2C=CC3)C)[H].C(=C\\C(O)=O)\\C(O)=O '
               'NAME: Methysergide maleate REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C(O[C@H]1[C@](O)([C@@H]2[C@H]3O[C@H]3CCN2[C@@H]1/C=C(/C=C/C4CC4)\\C)C)CC(C)C '
               'NAME: JBIR-102 REASON: WRONGLY CLASSIFIED Matches essential '
               'fatty acid criteria with multiple cis double bonds and a '
               'carboxylic acid group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C[C@H]1O[C@@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O '
               'NAME: (8R)-hydroxy-(11R,12R)-epoxyicosa-(5Z,9E,14Z)-trienoic '
               'acid REASON: WRONGLY CLASSIFIED Matches essential fatty acid '
               'criteria with multiple cis double bonds and a carboxylic acid '
               'group\n'
               ' * SMILES: '
               'OC12C(C3C(O)(CC(C3)(C)C)C=C1CO)(CC2OC(=O)C=4C(=CC(O)=CC4O)C)C '
               'NAME: 13-Hydroxydihydromelleolide REASON: WRONGLY CLASSIFIED '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               '[H+].[H+].[H+].[H+].[O-]C(=O)\\C=C/C([O-])=O.[O-]C(=O)\\C=C/C([O-])=O.CN1CCN(CC1)c1nc2cc(ccc2n2cccc12)C(F)(F)F '
               'NAME: '
               '4-(4-methylpiperazin-1-yl)-7-(trifluoromethyl)pyrrolo[1,2-a]quinoxaline '
               'dimaleate REASON: WRONGLY CLASSIFIED Matches essential fatty '
               'acid criteria with multiple cis double bonds and a carboxylic '
               'acid group\n'
               ' * SMILES: '
               'C1[C@@H]([C@H]([C@H]([C@@H]1O)/C=C/[C@H](CCCCC)O)C/C=C\\CCCC(O)=O)O '
               'NAME: 8-epi-prostaglandin F2alpha REASON: WRONGLY CLASSIFIED '
               'Matches essential fatty acid criteria with multiple cis double '
               'bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C1[C@@]2(C(=O)O[C@]34[C@@H]2O[C@@H]([C@@H]1CO)[C@@H](C[C@@H](CCC[C@H]3C=C(C(=O)O)[C@H](C4)C)C)C)C '
               'NAME: Okilactomycin A REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'C[C@H](CCC[C@H](C)C(O)=O)[C@H]1CC[C@H]2C3=CC[C@H]4CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C '
               'NAME: Delta(1),Delta(7)-dafachronic acid REASON: WRONGLY '
               'CLASSIFIED Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C1C(OC(=O)[C@@]2(C[C@H]3[C@@]4(C(=C(O)C([C@]4(O)C[C@H]3C2)=O)C)C)C)=C5C[C@H]6C[C@@](C[C@H]6[C@@]5(C1=C)C)(C(=O)O)C '
               'NAME: Sterhirsutin E REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'C[C@]12[C@H]3[C@H](C(O)=O)C4=C[C@H]5CC[C@@]4([C@@H](O)C5=C)[C@@]3(CC[C@@H]1O)OC2=O '
               'NAME: Antheridiogen-An REASON: WRONGLY CLASSIFIED Matches '
               'essential fatty acid criteria with multiple cis double bonds '
               'and a carboxylic acid group\n'
               ' * SMILES: '
               'O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)OC)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)C '
               'NAME: methyl 20(21)-dehydrolucidenate A REASON: WRONGLY '
               'CLASSIFIED Matches essential fatty acid criteria with multiple '
               'cis double bonds and a carboxylic acid group\n'
               'False negatives: SMILES: '
               'CCCCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC '
               'NAME: '
               '1-eicosanoyl-2-[(7Z,10Z,13Z,16Z,19Z)-docosapentaenoyl]-sn-glycero-3-phosphocholine '
               'REASON: MISSED Carbon chain length 50 is outside typical range '
               'for essential fatty acids (18 to 36 carbons)\n'
               ' * SMILES: CCCCC\\C=C/CCC(O)=O NAME: cis-4-decenoic acid '
               'REASON: MISSED Insufficient number of cis double bonds, found '
               '1\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/CCCCCC(O)=O NAME: '
               'all-cis-7,10,13-hexadecatrienoic acid REASON: MISSED Carbon '
               'chain length 16 is outside typical range for essential fatty '
               'acids (18 to 36 carbons)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-octadecanoyl-2-[(6Z,9Z,12Z)-octadecatrienoyl]-sn-glycero-3-phosphocholine '
               'REASON: MISSED Carbon chain length 44 is outside typical range '
               'for essential fatty acids (18 to 36 carbons)\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCCCCCC(O)=O '
               'NAME: (23Z,26Z,29Z,32Z)-octatriacontatetraenoic acid REASON: '
               'MISSED Carbon chain length 38 is outside typical range for '
               'essential fatty acids (18 to 36 carbons)\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCO[C@H](COC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC)COP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: 1-arachidonyl-2-palmityl-sn-glycero-3-phosphocholine '
               'REASON: MISSED Carbon chain length 44 is outside typical range '
               'for essential fatty acids (18 to 36 carbons)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Missing terminal carboxylic '
                                               'acid group (C(=O)O)'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Missing terminal carboxylic '
                                               'acid group (C(=O)O)'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Missing terminal carboxylic '
                                               'acid group (C(=O)O)'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Insufficient number of cis '
                                               'double bonds, found 0'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Missing terminal carboxylic '
                                               'acid group (C(=O)O)'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Insufficient number of cis '
                                               'double bonds, found 0'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Insufficient number of cis '
                                               'double bonds, found 0'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Missing terminal carboxylic '
                                               'acid group (C(=O)O)'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'Insufficient number of cis '
                                               'double bonds, found 0'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Insufficient number of cis '
                                               'double bonds, found 0'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCC\\C=C/CCC(O)=O',
                                      'name': 'cis-4-decenoic acid',
                                      'reason': 'Insufficient number of cis '
                                                'double bonds, found 1'}],
    'attempt': 1,
    'reasoning': 'The previous attempt to classify essential fatty acids was '
                 'promising but faced certain challenges with false negatives '
                 'and false positives. The issues can be addressed by refining '
                 'the criteria and examples used in the function.\n'
                 '\n'
                 '### Issues Analysis:\n'
                 '\n'
                 '1. **Carbon Chain Length**: The previous definition '
                 'considered essential fatty acids to have 18 to 36 carbon '
                 'atoms. However, examples in the list exceed this range '
                 '(e.g., with chain lengths of 38 to 50 carbons), and some are '
                 'shorter (e.g., 16 carbons).\n'
                 '\n'
                 '2. **Cis Double Bond Count**: While the function requires '
                 "multiple cis double bonds, it didn't account for fatty acids "
                 'with fewer or different configurations of double bonds being '
                 'valid. Some classifications rely heavily on other criteria '
                 'instead of just the count of cis double bonds, such as '
                 'specific functional groups (e.g., phosphocholine).\n'
                 '\n'
                 '3. **False Positives**: Some molecules were misclassified as '
                 'they have multiple double bonds and carboxylic acids but '
                 'belong to a different class structurally (e.g., lactones, '
                 'aromatics).\n'
                 '\n'
                 '4. **Use of Aromatic Bonds in SMARTS**: Check to exclude '
                 'aromatic rings which were previously misclassified due to '
                 'double bonds.\n'
                 '\n'
                 '### Proposed Improvements:\n'
                 '\n'
                 '1. **Adjust Carbon Chain Length Criteria**: Broaden '
                 'allowable chain lengths, or remove this constraint if it '
                 'does not help differentiate true positives and false '
                 'positives.\n'
                 '\n'
                 '2. **Include Specific Patterns for Other Functional '
                 'Groups**: In addition to carboxylic acids, examine other '
                 'functional attachments like phosphocholines, which are '
                 'present in the example list.\n'
                 '\n'
                 '3. **Exclude Aromatic Double Bonds**: Specifically exclude '
                 'aromatic molecules from the match to avoid '
                 'misclassification.\n'
                 '\n'
                 '4. **Review the SMARTS Patterns**: Adjust to specifically '
                 'target essential fatty acids and exclude non-fatty acid '
                 'structures more effectively.\n'
                 '\n'
                 'Here is the updated code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 41,
    'num_false_positives': 100,
    'num_true_negatives': 736,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.2907801418439716,
    'recall': 0.9761904761904762,
    'f1': 0.448087431693989,
    'accuracy': 0.8849658314350797,
    'negative_predictive_value': 0.9986431478968792}