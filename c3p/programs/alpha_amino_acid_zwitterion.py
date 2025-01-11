"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Zwitterion pattern for alpha-amino acids:
    # Positive nitrogen attached to one alpha carbon, the same carbon attached to a carboxylate
    # Pattern for a positively charged nitrogen attached to an alpha carbon
    pos_nitrogen_pattern = Chem.MolFromSmarts("[NH1+,NH2+,NH3+]")  # Consideration for different positive nitrogen

    # Pattern for carboxylate attached to the alpha carbon
    carboxylate_pattern = Chem.MolFromSmarts("[C](=O)[O-]")

    pos_nitrogen_matches = mol.HasSubstructMatch(pos_nitrogen_pattern)
    carboxylate_matches = mol.HasSubstructMatch(carboxylate_pattern)

    if pos_nitrogen_matches and carboxylate_matches:
        return True, "Contains alpha carbon with NH3+ and COO-, indicating alpha-amino-acid zwitterion"
    else:
        return False, "No zwitterionic alpha carbon with NH3+ and COO- found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78608',
                          'name': 'alpha-amino-acid zwitterion',
                          'definition': 'An amino acid-zwitterion obtained by '
                                        'transfer of a proton from the carboxy '
                                        'to the amino group of any alpha-amino '
                                        'acid; major species at pH 7.3.',
                          'parents': ['CHEBI:35238'],
                          'xrefs': ['MetaCyc:Alpha-Amino-Acids'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5725490196078432 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C1(=CNC2=C1C=CC(=C2)Cl)C[C@@H](C(=O)[O-])[NH3+] NAME: '
               '6-chloro-L-tryptophan zwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: C1(=CNC2=C1C=CC(=C2Cl)Cl)C[C@@H](C(=O)[O-])[NH3+] '
               'NAME: 6,7-dichloro-L-tryptophan zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: BrC1=CC(C[C@@H](C([O-])=O)[NH3+])=CC=C1O NAME: '
               '3-bromo-L-tyrosine zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: [NH3+][C@@H](C[C@@H](O)C=O)C([O-])=O NAME: '
               'L-4-hydroxyglutamate semialdehyde zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+][C@@H](Cc1c[nH]c2ccccc12)C([O-])=O NAME: '
               'L-tryptophan zwitterion REASON: CORRECT Contains alpha carbon '
               'with NH3+ and COO-, indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: [O-]C([C@H](CCCCN)[NH3+])=O NAME: L-lysine '
               'zwitterion REASON: CORRECT Contains alpha carbon with NH3+ and '
               'COO-, indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])CCCCCCSC NAME: '
               'L-tetrahomomethionine zwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: CC(=O)[C@H]([NH3+])C([O-])=O NAME: '
               'L-2-amino-3-oxobutanoic acid zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'O1[C@@H]([C@H]([C@H]([C@@H]([C@H]1CO)O)O)O)C2=C(C3=C(N2)C=CC=C3)C[C@@H](C(=O)[O-])[NH3+] '
               "NAME: 2'-alpha-mannosyl-L-tryptophan zwitterion REASON: "
               'CORRECT Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: C(\\S(C[C@@H](C([O-])=O)[NH3+])=O)=C/C NAME: '
               'S-1-propenyl-L-cysteine sulfoxide zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'C(=O)([O-])[C@@H]([NH3+])CSC/C=C(/CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)\\C '
               'NAME: S-[(2E,6E,10E)-geranylgeranyl]-L-cysteine zwitterion '
               'REASON: CORRECT Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               '[C@H]([C@H](C[C@H](/C=C/CCCC/C=C/CCCCCCC)O)O)(C(=O)[O-])[NH3+] '
               'NAME: presphingofungin zwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: '
               'N1(C([C@H](C1)NC(/C(/C=2C=CC(OCC[C@@H](C(=O)[O-])[NH3+])=CC2)=N\\O)=O)=O)[C@@H](C(O)=O)C3=CC=C(C=C3)O '
               'NAME: isonocardicin A zwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: CC[C@H]([NH3+])C([O-])=O NAME: '
               'L-alpha-aminobutyrate zwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: '
               '[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)[O-])=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@@H](C([O-])=O)[NH3+])[NH3+])=O '
               'NAME: isonocardicin C dizwitterion REASON: CORRECT Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])CC=1C=C(C(=CC1)O)OC NAME: '
               '3-O-methyldopa zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: '
               'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H]([NH3+])C([O-])=O)C([O-])=O '
               'NAME: hercynylcysteine sulfoxide zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+][C@@H](CC1=CC=C(O)C(Cl)=C1)C([O-])=O NAME: '
               '3-chloro-L-tyrosine zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])CCCCCCCCSC NAME: '
               'L-hexahomomethionine zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: C(=O)([C@@H]([NH3+])CCS(C)(=O)=O)[O-] NAME: '
               'methionine sulfone zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: CC(=O)CC([NH3+])C([O-])=O NAME: '
               '2-amino-4-oxopentanoic acid zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: C1(=CNC2=C1C=C(C=C2)Cl)C[C@@H](C(=O)[O-])[NH3+] '
               'NAME: 5-chloro-L-tryptophan zwitterion REASON: CORRECT '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: NC(=O)NOC[C@H]([NH3+])C([O-])=O NAME: '
               'O-ureido-L-serine zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])CCCCCSC NAME: '
               'L-trihomomethionine zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: [H]C(=O)C[C@H]([NH3+])C([O-])=O NAME: L-aspartic '
               'acid 4-semialdehyde zwitterion REASON: CORRECT Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               'False positives: SMILES: '
               'CCCCCCCCCCCCCCCC\\C=C/OC[C@H](COP([O-])(=O)OC[C@H]([NH3+])C([O-])=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-(1Z-octadecenyl)-2-arachidonoyl-sn-glycero-3-phosphoserine(1-) '
               'REASON: WRONGLY CLASSIFIED Contains alpha carbon with NH3+ and '
               'COO-, indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+]CCCCNC(=O)CC[C@H]([NH3+])C([O-])=O NAME: '
               'gamma-L-glutamylputrescinium(1+) REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: CSCCCCCCC([NH3+])C([O-])=O NAME: '
               'tetrahomomethionine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               '[NH3+][C@H](CCOc1ccc(cc1)C(=N\\[O-])\\C(=O)N[C@H]1CN([C@@H](C([O-])=O)c2ccc(O)cc2)C1=O)C([O-])=O '
               'NAME: nocardicin A(2-) REASON: WRONGLY CLASSIFIED Contains '
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: Nc1ccc(O)cc1C(=O)C[C@H]([NH3+])C([O-])=O NAME: '
               '5-hydroxy-L-kynurenine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'IC1=CC(OC2=CC=C(C[C@H]([NH3+])C([O-])=O)C=C2)=CC(I)=C1O NAME: '
               "3',5'-diiodo-L-thyronine REASON: WRONGLY CLASSIFIED Contains "
               'alpha carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: C([C@@H](CCP(=O)(C)[O-])[NH3+])([O-])=O NAME: '
               '(2R)-glufosinate zwitterion(1-) REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'C(=O)([C@@H](NC(=O)CC[C@@H](C(=O)[O-])[NH3+])[C@H](O)C)[O-] '
               'NAME: L-gamma-Glu-L-Thr(1-) REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               '[NH3+]C[C@@H](CC[C@H]([NH3+])C([O-])=O)OP([O-])([O-])=O NAME: '
               'erythro-5-phosphonatoooxy-L-lysinium(1-) REASON: WRONGLY '
               'CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: CSCCCCCCCC([NH3+])C([O-])=O NAME: '
               'pentahomomethionine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'C=1(C[C@@H](C(=O)[O-])[NH3+])C=C(C(OC2=CC(I)=C(O)C=C2)=CC1)I '
               "NAME: 3,3'-diiodo-L-thyronine zwitterion REASON: WRONGLY "
               'CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               '[NH3+][C@@H](CCCNC(=[NH2+])NP([O-])([O-])=O)C([O-])=O NAME: '
               'N(omega)-phosphonato-L-arginine REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCCCC)(OC[C@H]([NH3+])C([O-])=O)([O-])=O '
               'NAME: 24:0-18:2-PS REASON: WRONGLY CLASSIFIED Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: '
               '[NH3+][C@@H](CCCCNN[C@@H](CCC([O-])=O)C([O-])=O)C([O-])=O '
               'NAME: 1-L-glutamo-2-N(6-)L-lysinohydrazine(2-) REASON: WRONGLY '
               'CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])CC=1C=CC(=CC1)OS([O-])(=O)=O '
               "NAME: O(4')-sulfo-L-tyrosinate(1-) REASON: WRONGLY CLASSIFIED "
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP([O-])(=O)OC[C@H]([NH3+])C([O-])=O '
               'NAME: 1-stearoyl-sn-glycero-3-phosphoserine(1-) REASON: '
               'WRONGLY CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: CNC(=O)CC[C@H]([NH3+])C([O-])=O NAME: '
               'N(5)-methyl-L-glutamine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+][C@@H](C[SeH])C([O-])=O NAME: '
               'L-selenocysteine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCCCCCC)(OC[C@H]([NH3+])C([O-])=O)([O-])=O '
               'NAME: 20:0-18:3-PS REASON: WRONGLY CLASSIFIED Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: NC(=O)NCCCC[C@H]([NH3+])C([O-])=O NAME: '
               'L-homocitrulline zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+][C@H](C([O-])=O)C(=O)COP([O-])([O-])=O NAME: '
               '(2S)-2-ammonio-3-oxo-4-(phosphonatooxy)butanoate(2-) REASON: '
               'WRONGLY CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: [NH3+]CCCC[C@H]([NH3+])C([O-])=O NAME: '
               'L-lysinium(1+) REASON: WRONGLY CLASSIFIED Contains alpha '
               'carbon with NH3+ and COO-, indicating alpha-amino-acid '
               'zwitterion\n'
               ' * SMILES: CCCC(SC[C@@H](C([O-])=O)[NH3+])(C)CCO NAME: '
               'S-(1-hydroxy-3-methylhexan-3-yl)-L-cysteine zwitterion REASON: '
               'WRONGLY CLASSIFIED Contains alpha carbon with NH3+ and COO-, '
               'indicating alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               '[NH3+][C@@H](CCC(=O)N[C@@H](CSS)C(=O)NCC([O-])=O)C([O-])=O '
               'NAME: S-sulfanylglutathionate(1-) REASON: WRONGLY CLASSIFIED '
               'Contains alpha carbon with NH3+ and COO-, indicating '
               'alpha-amino-acid zwitterion\n'
               ' * SMILES: '
               'C[NH+](C)C(CCc1ncc(C[C@H]([NH3+])C([O-])=O)[nH]1)C([O-])=O '
               'NAME: 2-[3-carboxylato-3-(dimethylammonio)propyl]-L-histidine '
               'dizwitterion REASON: WRONGLY CLASSIFIED Contains alpha carbon '
               'with NH3+ and COO-, indicating alpha-amino-acid zwitterion\n'
               'False negatives: SMILES: '
               '[Se]([C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(C)=O)(C2=NC(=CN2)C[C@@H](C([O-])=O)[N+](C)(C)C)=O '
               'NAME: hercynyl-SeGlcNAc selenoxide REASON: MISSED No '
               'zwitterionic alpha carbon with NH3+ and COO- found\n'
               ' * SMILES: [O-]C(=O)C([NH2+]C)(C)C NAME: '
               'N-methyl-aminoisobutyrate zwitterion REASON: MISSED No '
               'zwitterionic alpha carbon with NH3+ and COO- found\n'
               ' * SMILES: O=C([O-])[C@@H](N)[C@@H](C)[NH3+] NAME: '
               '(2S,3R)-diazaniumylbutanoate zwitterion REASON: MISSED No '
               'zwitterionic alpha carbon with NH3+ and COO- found\n'
               ' * SMILES: O=C([O-])[C@@H](N)[C@H](C)[NH3+] NAME: '
               '(2S,3S)-diazaniumylbutanoate zwitterion REASON: MISSED No '
               'zwitterionic alpha carbon with NH3+ and COO- found\n'
               ' * SMILES: C[NH+](C)[C@@H](Cc1c(S)ncn1C)C([O-])=O NAME: '
               'ovothiol C zwitterion REASON: MISSED No zwitterionic alpha '
               'carbon with NH3+ and COO- found\n'
               ' * SMILES: C(C1[NH2+]CSC1)(=O)[O-] NAME: thioproline '
               'zwitterion REASON: MISSED No zwitterionic alpha carbon with '
               'NH3+ and COO- found\n'
               ' * SMILES: C[NH2+][C@@H](Cc1c(S)ncn1C)C([O-])=O NAME: ovothiol '
               'B zwitterion REASON: MISSED No zwitterionic alpha carbon with '
               'NH3+ and COO- found\n'
               ' * SMILES: O[C@H]1C[NH2+][C@@H](C1)C([O-])=O NAME: '
               'trans-4-hydroxy-L-proline zwitterion REASON: MISSED No '
               'zwitterionic alpha carbon with NH3+ and COO- found\n'
               ' * SMILES: [C@H]1(CCC[NH+]1C)C([O-])=O NAME: N-methylproline '
               'zwitterion REASON: MISSED No zwitterionic alpha carbon with '
               'NH3+ and COO- found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No zwitterionic alpha carbon '
                                               'with NH3+ and COO- found'}],
    'sample_false_negatives': [   {   'smiles': '[Se]([C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(C)=O)(C2=NC(=CN2)C[C@@H](C([O-])=O)[N+](C)(C)C)=O',
                                      'name': 'hercynyl-SeGlcNAc selenoxide',
                                      'reason': 'No zwitterionic alpha carbon '
                                                'with NH3+ and COO- found'}],
    'attempt': 3,
    'reasoning': 'The previous program defined a SMARTS pattern to identify an '
                 'alpha-amino acid zwitterion as having an alpha carbon '
                 'connected to both an ammonium group (`[NH3+]`) and a '
                 'carboxylate group (`[C(=O)[O-]]`). However, the outcomes '
                 'highlighted several issues that led to false positives and '
                 'false negatives.\n'
                 '\n'
                 '**Issues and Improvements:**\n'
                 '\n'
                 '1. **Overly Specific Pattern:**\n'
                 '   - The SMARTS pattern is too specific and does not account '
                 'for cases where the structure may vary, while still being an '
                 'alpha-amino acid zwitterion.\n'
                 '   - For instance, it may not identify molecules where the '
                 'ammonium or carboxylate groups are slightly altered or when '
                 'connected in a slightly different configuration.\n'
                 '\n'
                 '2. **Misclassification of Complex Structures:**\n'
                 '   - Complex molecules with additional functional groups or '
                 'differing backbone atoms might still classify as zwitterions '
                 'but were not recognized by the simple pattern.\n'
                 '\n'
                 '3. **Improved Pattern Matching:**\n'
                 '   - We need a more flexible SMARTS pattern to account for '
                 'variations in molecular structure while conserving the '
                 'zwitterionic characteristic.\n'
                 '   - We could split pattern checking into distinct checks '
                 'for the presence of core zwitterionic features rather than a '
                 'single pattern:\n'
                 '     - Presence of a positively charged nitrogen atom `[N+]` '
                 'attached to an alpha carbon.\n'
                 '     - Presence of a negatively charged carboxylate '
                 '`C(=O)[O-]` group also connected to the alpha carbon.\n'
                 '\n'
                 '4. **Avoid False Positives:**\n'
                 '   - Filter out structures that have additional charged '
                 'groups that classify them as something else entirely, or '
                 'have more complex structural elements like phosphates that '
                 'are not typical of simple zwitterions.\n'
                 '\n'
                 "Here's an improved version of the program using the detailed "
                 'observations:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 81,
    'num_false_positives': 100,
    'num_true_negatives': 16380,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.44751381215469616,
    'recall': 0.9878048780487805,
    'f1': 0.6159695817490495,
    'accuracy': 0.9939017026929114,
    'negative_predictive_value': 0.9999389536658324}