"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded patterns for detecting common bioconjugate motifs
    patterns = {
        "peptide_bond": Chem.MolFromSmarts("N[C@H](C)C(=O)"),
        "disulfide_bond": Chem.MolFromSmarts("S-S"),
        "glutathione_like": Chem.MolFromSmarts("N[C@H](CC(=O)NCCS)C(=O)NCC(=O)O"),
        "coenzyme_a_linkage": Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)"),
        "thioester_bond": Chem.MolFromSmarts("C(=O)S"),
        "ester_bond": Chem.MolFromSmarts("C(=O)O"),
        "amide_bond": Chem.MolFromSmarts("C(=O)N"),
        "heterocyclic_nitrogen": Chem.MolFromSmarts("n")
    }
    
    # Track which patterns we match
    matched_patterns = set()
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.add(name)

    # Determine if there are at least two distinctive substructures suggesting bioconjugation
    if len(matched_patterns) >= 2:
        return True, f"Contains patterns: {', '.join(matched_patterns)}"

    # For failed cases, provide reasoning based on detected patterns
    if matched_patterns:
        return False, f"Partially matched patterns, only found: {', '.join(matched_patterns)}"

    # Provide explanation when no patterns are matched
    return False, "No definitive bioconjugate patterns found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64985',
                          'name': 'bioconjugate',
                          'definition': 'A molecular entity consisting of at '
                                        'least 2 biological molecules '
                                        'covalently linked together.',
                          'parents': ['CHEBI:50860'],
                          'xrefs': ['Wikipedia:Bioconjugate'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.029411764705882353 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: S(SS(O)(=O)=O)C[C@H](N)C(O)=O NAME: '
               '3-(Sulfodisulfanyl)-L-alanine REASON: CORRECT Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S(SC[C@H](N)C(O)=O)C[C@H](C(=O)N1[C@@H](CCC1)C(O)=O)C NAME: '
               'Captopril-cysteine disulfide REASON: CORRECT Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               'False positives: SMILES: '
               'S1SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(NC(=O)[C@@H](NC(=O)CC12CCCCC2)CC3=CC=C(OC)C=C3)[C@H](CC)C)CCC(=O)N)CC(=O)N)C(=O)N4[C@@H](CCC4)C(=O)N[C@@H](CCCN)C(=O)NCC(=O)N '
               'NAME: Mpomeovt REASON: WRONGLY CLASSIFIED Contains patterns: '
               'peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'CC(CSSC[C@@H](C)C(=O)N1CCC[C@H]1C(O)=O)C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: captopril disulfide REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SCC(NC(=O)C2NC(=O)C(NC(=O)C(CC3=CC=C(O)C=C3)NC(C(NC(C4(CC(NC(C(C1)NC(=O)C(NC(=O)/C(/NC(=O)C(NC(=O)C(NC(=O)/C(/NC(=O)/C(/NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(N)CC(=O)O)C)C(CC)C)=C/C)=C/C)C(CC)C)C(CC)C)=C/C)C(CC)C)=O)C(=O)NC(C(=O)NC(C(N4)=O)CC(C)C)C(C)C)CSC2)=O)C(C)C)=O)CCCCN)C(=O)O '
               'NAME: Paenithopeptin C REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'C(N[C@H]1CSSC[C@H](NC([C@H](C(C)C)NC([C@H](CCCCN)NC([C@H](NC([C@@H](NC1=O)CC2=CC=C(O)C=C2)=O)CC3=CNC4=C3C=CC=C4)=O)=O)=O)C(N[C@@H]([C@@H](C)O)C(N)=O)=O)(=O)[C@](N)([H])CC=5C=C6C=CC=CC6=CC5 '
               'NAME: lanreotide REASON: WRONGLY CLASSIFIED Contains patterns: '
               'peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S(SCCNC(=O)[C@@H](NC(=O)CN)C)CCNC(=O)[C@@H](NC(=O)CN)C.FC(F)(F)C(O)=O.FC(F)(F)C(O)=O '
               'NAME: L-Alaninamide, glycyl-N-(2-mercaptoethyl)-, '
               "2,2'-disulfide bis(trifluoroacetate) REASON: WRONGLY "
               'CLASSIFIED Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: C(C(CCSSCCC(C(=O)O)N)N)(O)=O NAME: homocystine '
               'REASON: WRONGLY CLASSIFIED Contains patterns: peptide_bond, '
               'disulfide_bond\n'
               ' * SMILES: C([C@H](N)C(=O)O)SS NAME: 3-disulfanyl-L-alanine '
               'REASON: WRONGLY CLASSIFIED Contains patterns: peptide_bond, '
               'disulfide_bond\n'
               ' * SMILES: '
               '[H][C@]1(CSSC[C@]2([H])NC(=O)[C@H](Cc3ccc(O)cc3)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@@H]3CSSC[C@H](NC(=O)[C@H](C)NC(=O)[C@@H]4CCCN4C(=O)[C@H](CC(N)=O)NC(=O)[C@H](CSSC[C@H](N)C(=O)N3)NC2=O)C(=O)N[C@@]([H])([C@@H](C)O)C(=O)NCC(=O)N1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O '
               'NAME: linaclotide REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: N[C@@H](CSSC[C@H](N)C([O-])=O)C([O-])=O NAME: '
               'L-cystine anion REASON: WRONGLY CLASSIFIED Contains patterns: '
               'peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'C[C@@H](O)[C@@H]1NC(=O)[C@H](CCCCN)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CSSC[C@H](NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc2ccccc2)NC1=O)[C@@H](C)O)C(O)=O)NC(=O)CNC(=O)[C@H](C)N '
               'NAME: somatostatin REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'N[C@@H](CCC(=O)N[C@H]1CSSC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(=O)NCCCNCCCCNC(=O)CNC1=O)C(O)=O '
               'NAME: trypanothione disulfide REASON: WRONGLY CLASSIFIED '
               'Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SC[C@@H]2N(C(=O)CNC(=O)[C@H](NC(=O)C3=NC=4C(=CC=CC4)C=C3)CSC([C@@H](N(C([C@H](C1)N(C(=O)CNC(=O)[C@H](NC(=O)C5=NC=6C(=CC=CC6)C=C5O)CSC(C(N(C2=O)C)CSC)=O)C)=O)C)CSC)=O)C '
               "NAME: 22'-Deoxythiocoraline REASON: WRONGLY CLASSIFIED "
               'Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: [NH3+][C@@H](CSS)C([O-])=O NAME: '
               '3-disulfanyl-L-alanine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S(=O)(C[C@@H]1N(C(=O)[C@H]2N(C(=O)CNC(=O)[C@H](NC(=O)C3=NC=4C(=CC=CC4)C=C3O)CSC(=O)[C@@H](N(C([C@H](CSSC2)N(C)C(CNC([C@@H](CSC1=O)NC(=O)C5=NC=6C(=CC=CC6)C=C5O)=O)=O)=O)C)CSC)C)C)C '
               'NAME: 12′-sulfoxythiocoraline REASON: WRONGLY CLASSIFIED '
               'Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: N[C@@H](CCCCNC(=O)CCCCC1CCSS1)C(O)=O NAME: '
               'N(6)-lipoyl-L-lysine REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               '[H]N[C@H]1CSSC[C@H](NC([C@@H](NC([C@@H](NC([C@@H](NC([C@@H](NC1=O)CC2=CC=C(O)C=C2)=O)CC3=CC=CC=C3)=O)CCC(N)=O)=O)CC(N)=O)=O)C(N4CCC[C@H]4C(N[C@@H](CCCN)C(NCC(N)=O)=O)=O)=O '
               'NAME: ornipressin REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SC[C@H](N(C(=O)CNC([C@@H](C1)NC(=O)C2=NC=3C(=CC=CC3)C=C2O)=O)C)C(=O)N([C@H](C(=O)OC)CSC)C '
               'NAME: Thiochondrilline C REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'CC1CC11N(C)C(=O)C2CSSCC(N(C)C(=O)C(C)NC(=O)C(COC1=O)NC(=O)c1nc3ccccc3cc1O)C(=O)N(C)C1(CC1C)C(=O)OCC(NC(=O)c1nc3ccccc3cc1O)C(=O)NC(C)C(=O)N2C '
               'NAME: SW 163C REASON: WRONGLY CLASSIFIED Contains patterns: '
               'peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SCC(NC(=O)C2NC(=O)C(NC(=O)C(CC3=CC=C(O)C=C3)NC(C(NC(C4(CC(NC(C(C1)NC(=O)C(NC(=O)/C(/NC(=O)C(NC(=O)C(NC(=O)/C(/NC(=O)/C(/NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(N)C)C)CC(=O)O)C)C(CC)C)=C/C)=C/C)C(CC)C)C(CC)C)=C/C)C(CC)C)=O)C(=O)NC(C(=O)NC(C(N4)=O)CC(C)C)C(C)C)CSC2)=O)C(C)C)=O)CCCCN)C(=O)O '
               'NAME: Paenithopeptin E REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SC[C@H]2NC(=O)[C@@H](NC(=O)[C@H](NC([C@@H](NC([C@@H](C1)NC2=O)=O)C(C)C)=O)CC(C)C)[C@H](CC)C '
               'NAME: Malformin A1 REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SCC2N(C(=O)C(NC(=O)C(NC(=O)C=3C=NC4=CC=CC=C4C3)COC(C(N(C(C(C1)N(C(=O)C(NC(=O)C(NC(=O)C=5C=NC6=CC=CC=C6C5)COC(C(N(C2=O)C)C(C)C)=O)C)C)=O)C)C(C)C)=O)C)C '
               'NAME: QN-Triostin A REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S(SCCNC(=O)CC[C@H](N)C(O)=O)CCNC=1C(=O)C=2[C@H]([C@]3(OC)N(C[C@@H]4N[C@H]34)C2C(=O)C1C)COC(=O)N '
               'NAME: '
               '7-N-(2-((2-(gamma-L-Glutamylamino)ethyl)dithio)ethyl)mitomycin '
               'C REASON: WRONGLY CLASSIFIED Contains patterns: peptide_bond, '
               'disulfide_bond\n'
               ' * SMILES: '
               'N[C@@H](CCC(=O)N[C@@H](CSSC[C@H](N)C(O)=O)C(=O)NCC(O)=O)C(O)=O '
               'NAME: L-cysteine glutathione disulfide REASON: WRONGLY '
               'CLASSIFIED Contains patterns: peptide_bond, disulfide_bond\n'
               ' * SMILES: '
               'S1SC[C@@H]2NC(=O)CNC(=O)[C@@H](NC(=O)[C@H]3N(C(=O)[C@@H](NC(=O)CNC(C[C@H](NC([C@@H](NC([C@H]4N(C2=O)CCC4)=O)CO)=O)C(=O)N[C@H](C(=O)N5[C@H](C(=O)NCC(=O)N[C@H](C(N[C@H](C(N[C@H](C(N6[C@H](C(N[C@H](C(N[C@H](C(N[C@@H](C1)C(=O)O)=O)C)=O)CC=7C8=C(C=CC=C8)NC7)=O)CCC6)=O)[C@H](O)C)=O)CC(=O)N)=O)CC=9C%10=C(C=CC=C%10)NC9)CCC5)[C@H](CC)C)=O)CC(C)C)CCC3)CC=%11C%12=C(C=CC=C%12)NC%11 '
               'NAME: BI-32169 REASON: WRONGLY CLASSIFIED Contains patterns: '
               'peptide_bond, disulfide_bond\n'
               ' * SMILES: S(SC(C(N)C(O)=O)(C)C)C(C(N)C(O)=O)(C)C NAME: '
               'Penicillamine disulfide REASON: WRONGLY CLASSIFIED Contains '
               'patterns: peptide_bond, disulfide_bond\n'
               'False negatives: SMILES: '
               'N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(O)=O NAME: L-Cys(IAN) '
               'REASON: MISSED Partially matched patterns, only found: '
               'peptide_bond\n'
               ' * SMILES: '
               'N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O '
               'NAME: gammaGluCys(IAN)Gly REASON: MISSED Partially matched '
               'patterns, only found: peptide_bond\n'
               ' * SMILES: '
               'N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(O)=O)C(O)=O '
               'NAME: gammaGluCys(IAN) REASON: MISSED Partially matched '
               'patterns, only found: peptide_bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N '
               'NAME: (5Z,8Z,11Z,14Z,17Z)-icosapentaenoyl-CoA REASON: MISSED '
               'No definitive bioconjugate patterns found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCCCCCCCCCCCC)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: palmitoyl-CoA REASON: MISSED No definitive bioconjugate '
               'patterns found\n'
               ' * SMILES: '
               '[C@H](C(NCC(=O)O)=O)(NC(CC[C@@H](C(=O)O)N)=O)CSC(CC(=O)[H])C '
               'NAME: S-(4-oxobutan-2-yl)glutathione REASON: MISSED Partially '
               'matched patterns, only found: peptide_bond\n'
               ' * SMILES: '
               '[C@H]1([C@H]([C@H](CC1=O)SC[C@H](NC(CC[C@H](N)C(=O)O)=O)C(=O)NCC(=O)O)/C=C/[C@H](CCCCC)O)C/C=C\\CCCC(=O)O '
               'NAME: (S)-PGA2-S-glutathione conjugate REASON: MISSED '
               'Partially matched patterns, only found: peptide_bond\n'
               ' * SMILES: '
               'COc1cccc2[nH]cc(CSC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(O)=O)c12 '
               'NAME: gammaGluCys(4MeOI3M)Gly REASON: MISSED Partially matched '
               'patterns, only found: peptide_bond\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA REASON: MISSED No '
               'definitive bioconjugate patterns found\n'
               ' * SMILES: S(=O)(CCCCNC(SC[C@H](N)C(O)=O)=S)C NAME: '
               'sulforaphane-cysteine REASON: MISSED Partially matched '
               'patterns, only found: peptide_bond\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No definitive bioconjugate '
                                               'patterns found'},
                                 {   'smiles': 'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                     'name': '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-)',
                                     'reason': 'Partially matched patterns, '
                                               'only found: ester_bond'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'Partially matched patterns, '
                                               'only found: amide_bond'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No definitive bioconjugate '
                                               'patterns found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'Partially matched patterns, '
                                               'only found: amide_bond'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No definitive bioconjugate '
                                               'patterns found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No definitive bioconjugate '
                                               'patterns found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No definitive bioconjugate '
                                               'patterns found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'Partially matched patterns, '
                                               'only found: amide_bond'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'Partially matched patterns, '
                                               'only found: '
                                               'heterocyclic_nitrogen'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'The previous program attempted to classify bioconjugates '
                 'using a set of predefined substructure patterns, such as '
                 'peptide bonds and disulfide bonds, and required the presence '
                 'of at least two distinct biological unit patterns to '
                 'classify a compound as a bioconjugate. However, the '
                 "program's F1 score was very low, indicating that it often "
                 'misclassified compounds.\n'
                 '\n'
                 'To improve the classification, we need to understand the '
                 'common structural motifs in the given examples. '
                 'Bioconjugates typically involve connections like peptide '
                 'bonds, disulfide bonds, and other specific linkages such as '
                 'ester or thioester bonds between distinct biopolymers or '
                 'biological units.\n'
                 '\n'
                 'The issues with the previous implementation likely stem '
                 'from:\n'
                 '1. Over-reliance on a small set of patterns that might not '
                 'cover the diversity of bioconjugates.\n'
                 '2. Insufficient recognition of unique linkages that '
                 'characterize bioconjugates.\n'
                 '\n'
                 'To address these issues, we can:\n'
                 '- Expand the set of substructure patterns to include more '
                 'specific linkages or structures commonly found in '
                 'bioconjugates.\n'
                 '- Consider not just the presence of multiple patterns but '
                 'also the context that indicates distinct biological units '
                 '(e.g., large scaffold elements connected by specific '
                 'linkages).\n'
                 '\n'
                 "Here's an updated version of the program that includes "
                 'additional patterns and checks for more diverse linkages:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 34,
    'num_false_positives': 100,
    'num_true_negatives': 215,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2537313432835821,
    'recall': 1.0,
    'f1': 0.40476190476190477,
    'accuracy': 0.7134670487106017,
    'negative_predictive_value': 1.0}