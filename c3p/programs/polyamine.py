"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Invalid SMILES string

    # Define amino group patterns: -NH2, -NH-, -NR2, [NH3+], [NH2+], etc.
    amino_patterns = [
        Chem.MolFromSmarts("[NX3;H2]"),  # -NH2
        Chem.MolFromSmarts("[NX3;H1]"),  # -NH-
        Chem.MolFromSmarts("[NX3;H0]"),  # -NR2
        Chem.MolFromSmarts("[NH3+]"),    # [NH3+]
        Chem.MolFromSmarts("[NH2+]"),    # [NH2+]
        Chem.MolFromSmarts("[N+](=O)[O-]"),  # Nitro group
        Chem.MolFromSmarts("[N+]([O-])=O"),  # Nitro group (alternative)
    ]

    # Count amino groups
    amino_count = 0
    for pattern in amino_patterns:
        if pattern is not None:  # Ensure the pattern is valid
            amino_count += len(mol.GetSubstructMatches(pattern))

    # Check if there are at least two amino groups
    if amino_count >= 2:
        return True, f"Contains {amino_count} amino groups"
    else:
        return False, f"Contains only {amino_count} amino groups (need at least 2)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:88061',
                          'name': 'polyamine',
                          'definition': 'Any organic amino compound that '
                                        'contains two or more amino groups.',
                          'parents': ['CHEBI:50047'],
                          'xrefs': ['Wikipedia:Polyamine'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7055555555555556 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: NCCCCN NAME: putrescine REASON: '
               'CORRECT Contains 2 amino groups\n'
               ' * SMILES: Nc1cccc2c(N)cccc12 NAME: naphthalene-1,5-diamine '
               'REASON: CORRECT Contains 2 amino groups\n'
               " * SMILES: CCNCCNCC NAME: N,N'-diethylethylenediamine REASON: "
               'CORRECT Contains 2 amino groups\n'
               ' * SMILES: Nc1nc(Cl)nc(Nc2cccc(c2)S([O-])(=O)=O)n1 NAME: '
               '2-amino-6-chloro-4-(3-sulfonatoanilino)-1,3,5-triazine REASON: '
               'CORRECT Contains 2 amino groups\n'
               ' * SMILES: CC(C)(C)Nc1nc(N)nc(Cl)n1 NAME: '
               'terbutylazine-desethyl REASON: CORRECT Contains 2 amino '
               'groups\n'
               ' * SMILES: NCCCNCCSP(O)(O)=O NAME: amifostine REASON: CORRECT '
               'Contains 2 amino groups\n'
               ' * SMILES: NCCCCCCN NAME: hexane-1,6-diamine REASON: CORRECT '
               'Contains 2 amino groups\n'
               ' * SMILES: CCNc1nc(Cl)nc(NC(C)(C)C#N)n1 NAME: cyanazine '
               'REASON: CORRECT Contains 2 amino groups\n'
               ' * SMILES: NC[C@@H](O)C[C@@H](O)CNC[C@H](O)C[C@H](O)CN NAME: '
               'pavettamine REASON: CORRECT Contains 3 amino groups\n'
               ' * SMILES: CCNc1nc(N)nc(Cl)n1 NAME: deisopropylatrazine '
               'REASON: CORRECT Contains 2 amino groups\n'
               ' * SMILES: NCCCCNCCCN NAME: spermidine REASON: CORRECT '
               'Contains 3 amino groups\n'
               ' * SMILES: Nc1ncnc(N)n1 NAME: 1,3,5-triazine-2,4-diamine '
               'REASON: CORRECT Contains 2 amino groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](C[C@@H]([C@H]([C@@H]1O)O[C@H]2O[C@@H](C(C[C@H]2N)=O)CO)N)N '
               "NAME: 4'-oxolividamine REASON: CORRECT Contains 3 amino "
               'groups\n'
               ' * SMILES: NCCCN[13CH2][13CH2][13CH2][13CH2]N NAME: '
               'spermidine-(butyl-(13)C4) REASON: CORRECT Contains 3 amino '
               'groups\n'
               ' * SMILES: CC(=O)NCCCCN NAME: N-acetylputrescine REASON: '
               'CORRECT Contains 2 amino groups\n'
               ' * SMILES: '
               'C(CN(CC1=NC=CC=C1)CC2=NC=CC=C2)N(CC3=NC=CC=C3)CC4=NC=CC=C4 '
               "NAME: N,N,N',N'-tetrakis(2-pyridylmethyl)ethylenediamine "
               'REASON: CORRECT Contains 2 amino groups\n'
               ' * SMILES: NCCCCN(CCCN)CCCN NAME: N(4)-aminopropylspermidine '
               'REASON: CORRECT Contains 4 amino groups\n'
               ' * SMILES: NCCCNCCCCN(CCCN)CCCN NAME: N(4)-aminopropylspermine '
               'REASON: CORRECT Contains 5 amino groups\n'
               ' * SMILES: CNCCCCN NAME: N-methylputrescine REASON: CORRECT '
               'Contains 2 amino groups\n'
               ' * SMILES: N(c1ccccc1)c1ccc(Nc2ccccc2)cc1 NAME: '
               "N,N'-diphenyl-1,4-phenylenediamine REASON: CORRECT Contains 2 "
               'amino groups\n'
               ' * SMILES: NCCNCCNCCNCCN NAME: tetraethylenepentamine REASON: '
               'CORRECT Contains 5 amino groups\n'
               ' * SMILES: '
               'N[C@@H](CSc1cc(C[C@H](N)C(O)=O)c(SC[C@H](N)C(O)=O)c(O)c1O)C(O)=O '
               "NAME: 2,5-S,S'-dicysteinyldopa REASON: CORRECT Contains 3 "
               'amino groups\n'
               ' * SMILES: '
               'C=1C(=CC(C)=CC1C)OC[C@@H](NC2=NC(=NC(=N2)N)C(C)(C)F)C NAME: '
               '(S)-triaziflam REASON: CORRECT Contains 2 amino groups\n'
               ' * SMILES: COc1nc(NC(C)C)nc(NC(C)C)n1 NAME: prometon REASON: '
               'CORRECT Contains 2 amino groups\n'
               ' * SMILES: CCNc1nc(NC(C)CC)nc(OC)n1 NAME: '
               "N-(butan-2-yl)-N'-ethyl-6-methoxy-1,3,5-triazine-2,4-diamine "
               'REASON: CORRECT Contains 2 amino groups\n'
               'False positives: SMILES: '
               'C[C@@H]1CN([C@@H](COC2=C(C=C(C=C2)NS(=O)(=O)C3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C '
               'NAME: '
               'N-[(4R,7R,8S)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]benzenesulfonamide '
               'REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               ' * SMILES: '
               'O=C(N[C@@H](CO)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1NC=NC1)CCC(O)=O '
               'NAME: His-Glu-Ser REASON: WRONGLY CLASSIFIED Contains 3 amino '
               'groups\n'
               ' * SMILES: '
               'O=C1N(O)CCOCCNC(=O)CCC(=O)N(O)CCOCCNC(CCC(N(CCOCCNC(CC1)=O)O)=O)=O '
               'NAME: Desferrioxamine Et3 REASON: WRONGLY CLASSIFIED Contains '
               '6 amino groups\n'
               ' * SMILES: '
               'S(C[C@H](N)C(=O)N[C@H]1[C@@H](O)C(N2C3=NC=NC(=C3N=C2)N(C)C)OC1CO)C '
               'NAME: Cystocin REASON: WRONGLY CLASSIFIED Contains 3 amino '
               'groups\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(C(=CC=C2)NC(=O)C3CC3)O[C@@H]1CN(C)S(=O)(=O)C)[C@@H](C)CO '
               'NAME: '
               'N-[(2S,3R)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl(methylsulfonyl)amino]methyl]-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]cyclopropanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               ' * SMILES: '
               'CC1=C(C(=NO1)C)NC(=O)NCC[C@@H]2CC[C@H]([C@@H](O2)CO)NC(=O)C3=CC=C(C=C3)F '
               'NAME: '
               'N-[(2R,3R,6S)-6-[2-[[[(3,5-dimethyl-4-isoxazolyl)amino]-oxomethyl]amino]ethyl]-2-(hydroxymethyl)-3-oxanyl]-4-fluorobenzamide '
               'REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)[C@@H]1O)CO)[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(=O)NC '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5R,6S)-3-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-hydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxy-N-methyl-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               ' * SMILES: '
               'C[C@@H]1CCCCO[C@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)NS(=O)(=O)C)[C@@H](C)CO)C)CN(C)C(=O)NC3=CC=C(C=C3)C(F)(F)F '
               'NAME: '
               '1-[[(3R,9R,10R)-12-[(2S)-1-hydroxypropan-2-yl]-16-(methanesulfonamido)-3,10-dimethyl-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-9-yl]methyl]-1-methyl-3-[4-(trifluoromethyl)phenyl]urea '
               'REASON: WRONGLY CLASSIFIED Contains 4 amino groups\n'
               ' * SMILES: '
               'O1[C@@H](O[C@H]2[C@@H](O)[C@H](O[C@@H](O[C@H]3[C@@H](O)[C@@H](O)C(O[C@@H]3CO)O)[C@H]2O)CO[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]4NC(=O)C)CO)[C@H](NC(=O)C)[C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O[C@@H]8O[C@H]([C@@H](O)[C@@H](O)[C@@H]8O)C)CO)[C@H](O[C@@H]9O[C@H]([C@@H](O)[C@@H](O)[C@@H]9O)C)[C@H]1CO '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4S,5R,6S)-6-[(2R,3S,4R,5R,6R)-5-acetamido-6-[[(2R,3S,4S,5R,6S)-4-[(2S,3R,4R,5S,6R)-3-acetamido-4-[(2R,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-6-(hydroxymethyl)-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-[(2R,3S,4R,5R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-2-yl]methoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               ' * SMILES: '
               'C[C@H]1O[C@H](OP([O-])(=O)OP([O-])(=O)OC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(\\C)CC\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C)[C@H](NC(C)=O)[C@@H](O[C@H]2O[C@H](CO)[C@H](O[C@H]3O[C@H](CO)[C@H](O[C@H]4O[C@H](CO)[C@H](O[C@H]5O[C@H](CO)[C@H](O[C@H]6O[C@H](CO)[C@H](O)[C@H](O)[C@H]6NC(C)=O)[C@H](O)[C@H]5NC(C)=O)[C@H](O)[C@H]4NC(C)=O)[C@H](O)[C@H]3NC(C)=O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1NC(C)=O '
               'NAME: '
               '[alpha-D-GalNAc-(1->4)]4-alpha-D-GalNAc-(1->3)-alpha-D-diNAcBac-tritrans,heptacis-undecaprenyl '
               'diphosphate(2-) REASON: WRONGLY CLASSIFIED Contains 7 amino '
               'groups\n'
               ' * SMILES: '
               'CCN1C[C@H]([C@@H](CN(C(=O)C2=C(C=CC(=C2)NC(=O)NC(C)C)OC[C@H]1C)C)OC)C '
               'NAME: '
               '1-[(4R,7R,8S)-5-ethyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-propan-2-ylurea '
               'REASON: WRONGLY CLASSIFIED Contains 4 amino groups\n'
               ' * SMILES: CCNC(=S)N(CC1=CC=CO1)CC2=CC3=CC4=C(C=C3NC2=O)OCCO4 '
               'NAME: '
               '3-ethyl-1-(2-furanylmethyl)-1-[(7-oxo-3,6-dihydro-2H-[1,4]dioxino[2,3-g]quinolin-8-yl)methyl]thiourea '
               'REASON: WRONGLY CLASSIFIED Contains 2 amino groups\n'
               ' * SMILES: '
               'C1CC(C1)NC(=O)C[C@H]2CC[C@@H]([C@H](O2)CO)NS(=O)(=O)C3=CC=C(C=C3)F '
               'NAME: '
               'N-cyclobutyl-2-[(2R,5S,6S)-5-[(4-fluorophenyl)sulfonylamino]-6-(hydroxymethyl)-2-oxanyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains 2 amino groups\n'
               ' * SMILES: '
               'O[C@@H]([C@H](NC(=O)CNC(=O)[C@@H](N)CCC(O)=O)C(O)=O)C NAME: '
               'Glu-Gly-Thr REASON: WRONGLY CLASSIFIED Contains 3 amino '
               'groups\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C(=CC=C2)NC(=O)NC3=CC4=C(C=C3)OCO4)O[C@@H]1CNC)[C@@H](C)CO '
               'NAME: '
               '1-(1,3-benzodioxol-5-yl)-3-[(2S,3S)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-2-(methylaminomethyl)-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]urea '
               'REASON: WRONGLY CLASSIFIED Contains 4 amino groups\n'
               ' * SMILES: '
               'O=C1N[C@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@H](C(NCCCC[C@H]1NC(=O)[C@@H](NC(=O)CNC(=O)C[C@H](O)CCCCCCCCCCCC)CCCCN)=O)CCCCN)CCCCN)CCCCN)CCCCN '
               "NAME: PAX3' REASON: WRONGLY CLASSIFIED Contains 13 amino "
               'groups\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O)O[C@@H]([C@H]1O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@@H]7O[C@@H]([C@@H](O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8O)CO)[C@H](O)[C@H]7NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains 4 amino groups\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@H]1CN(C)C(=O)N[C@H](C)C(C)C)C4=CC=CC=C4N2C)[C@@H](C)CO '
               'NAME: LSM-30620 REASON: WRONGLY CLASSIFIED Contains 3 amino '
               'groups\n'
               ' * SMILES: '
               'CC1=CC=CC=C1OCC2=C(C=CC(=C2)C3NC4=CC=CC=C4C(=O)N3CC5=CC=CO5)OC '
               'NAME: '
               '3-(2-furanylmethyl)-2-[4-methoxy-3-[(2-methylphenoxy)methyl]phenyl]-1,2-dihydroquinazolin-4-one '
               'REASON: WRONGLY CLASSIFIED Contains 2 amino groups\n'
               ' * SMILES: '
               'OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C)CC(O)=O)[C@H](CC)C '
               'NAME: Ala-Asp-Ile REASON: WRONGLY CLASSIFIED Contains 3 amino '
               'groups\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@@H]8[C@H](O)[C@@H](O[C@@H]([C@H]8O)CO[C@H]9O[C@@H]([C@@H](O)[C@H](O)[C@@H]9O[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O)[C@H](O)[C@H]%11O)CO)[C@H](O)[C@H]%10NC(=O)C)CO)CO)O[C@H]%12[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]%12CO)O[C@H]%13[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%13CO)O '
               'NAME: CID 91856537 REASON: WRONGLY CLASSIFIED Contains 7 amino '
               'groups\n'
               ' * SMILES: '
               '[H][C@](C)(CC)c1nc2c(CCCNC(N)=N)[nH]c(cn2c1=O)-c1c[nH]c2ccccc12 '
               'NAME: Cypridina luciferin REASON: WRONGLY CLASSIFIED Contains '
               '2 amino groups\n'
               ' * SMILES: '
               'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CCCC3)O[C@@H]1CN(C)CC4=CC=CC=N4)[C@H](C)CO '
               'NAME: '
               '(2R)-2-[(4R,5S)-8-(1-cyclopentenyl)-4-methyl-5-[[methyl(2-pyridinylmethyl)amino]methyl]-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-2-yl]-1-propanol '
               'REASON: WRONGLY CLASSIFIED Contains 2 amino groups\n'
               ' * SMILES: '
               'C1=CC=C(C=C1)C(C2=CC=CC=C2)C(=O)NC(C(Cl)(Cl)Cl)NC(=S)NC3=CC(=C(C=C3)F)[N+](=O)[O-] '
               'NAME: '
               '2,2-diphenyl-N-[2,2,2-trichloro-1-[[(4-fluoro-3-nitroanilino)-sulfanylidenemethyl]amino]ethyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains 4 amino groups\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C(=CC=C2)NS(=O)(=O)C3=CC=C(C=C3)C)O[C@@H]1CN(C)C)[C@H](C)CO '
               'NAME: '
               'N-[(2S,3S)-2-[(dimethylamino)methyl]-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]-4-methylbenzenesulfonamide '
               'REASON: WRONGLY CLASSIFIED Contains 3 amino groups\n'
               'False negatives: SMILES: [NH3+]CCCCN(C(C)=O)O NAME: '
               'N(1)-acetyl-N(1)-hydroxyputrescine(1+) REASON: MISSED Contains '
               'only 1 amino groups (need at least 2)\n'
               ' * SMILES: C(C[NH3+])CC[NH2+]CCC([O-])=O NAME: putreanine(1+) '
               'REASON: MISSED Contains only 0 amino groups (need at least 2)\n'
               ' * SMILES: '
               '[H][C@@]12Cc3c[nH]c4cccc(c34)[C@@]1([H])C[C@@H](CSC)CN2CCC '
               'NAME: pergolide REASON: MISSED Contains only 1 amino groups '
               '(need at least 2)\n'
               ' * SMILES: [NH3+]CCCCNO NAME: N-hydroxyputrescine(1+) REASON: '
               'MISSED Contains only 1 amino groups (need at least 2)\n'
               ' * SMILES: [H][C@@]12Cc3c[nH]c4cccc(c34)[C@@]1([H])CCCN2 NAME: '
               'ergoline REASON: MISSED Contains only 1 amino groups (need at '
               'least 2)\n'
               ' * SMILES: [NH2+](CCCC[NH2+]CCCC[NH3+])CCC[NH3+] NAME: '
               'N(1)-aminopropylhomospermidine(+4) REASON: MISSED Contains '
               'only 0 amino groups (need at least 2)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1=CC=C2C(=C1)C3=C(C4=C(C=CC=N4)C=C3)NS2(=O)=O',
                                     'name': '5H-quinolino[8,7-c][1,2]benzothiazine '
                                             '6,6-dioxide',
                                     'reason': 'Contains only 1 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'O([C@@H]1C([C@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@]([C@H](O)C4)(CCC(C5)(C)C)C(O)=O)[H])=CC3)C)(CC2)C)[H])(CC1)C)[H])(C)C)[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO',
                                     'name': 'Echinocystic acid 3-glucoside',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'O=C1O[C@](C(=O)OC)(CC2=CC(=C(O)C=C2)CC(=O)C)C(=C1O)C3=CC=C(OC)C=C3',
                                     'name': 'Versicolactone A',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O',
                                     'name': 'PE(18:3(9Z,12Z,15Z)/18:1(11Z))',
                                     'reason': 'Contains only 1 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'COc1cc(O)c(CC=C(C)C)cc1C(=O)[C@@H](O)Cc1ccc(O)cc1',
                                     'name': 'lespeflorin C3',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'NC(=O)OP(O)(O)=O',
                                     'name': 'carbamoyl phosphate',
                                     'reason': 'Contains only 1 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'O=C1OC(=CC(=C1C)O)C[C@@H](O)CC/C=C/C',
                                     'name': '6-(2′S-hydroxy-5′E-ene-1′-heptyl)-4-hydroxy-3-methyl-2H-pyran-2-one',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': '[H][C@]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O[C@H]1CCC(C)(C)[C@]2([H])[C@H](OC(C)=O)[C@H](O)[C@@]3(C)O[C@](C)(CC(=O)[C@]3([H])[C@@]12C)C=C',
                                     'name': 'forskoditerpenoside C',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': 'O=C(OC[C@@H](O)CO)C1=C(OC)C(=C(C(=O)OC)C=C1C)C(=O)C2=C(O)C=C(O)C=C2O',
                                     'name': 'Wentiphenone A',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'},
                                 {   'smiles': '[H][C@@]1(O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)Oc1c[nH]c2ccc(Br)c(Cl)c12',
                                     'name': '5-bromo-4-chloro-3-indolyl '
                                             'alpha-D-galactoside',
                                     'reason': 'Contains only 0 amino groups '
                                               '(need at least 2)'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12Cc3c[nH]c4cccc(c34)[C@@]1([H])C[C@@H](CSC)CN2CCC',
                                      'name': 'pergolide',
                                      'reason': 'Contains only 1 amino groups '
                                                '(need at least 2)'},
                                  {   'smiles': '[H][C@@]12Cc3c[nH]c4cccc(c34)[C@@]1([H])CCCN2',
                                      'name': 'ergoline',
                                      'reason': 'Contains only 1 amino groups '
                                                '(need at least 2)'}],
    'attempt': 1,
    'reasoning': 'The previous code had a few issues that led to incorrect '
                 'classifications:\n'
                 '\n'
                 "1. **Charged Amino Groups**: The code didn't account for "
                 'charged amino groups (e.g., `[NH3+]`), which are common in '
                 'polyamines. This led to false negatives where molecules with '
                 'charged amino groups were incorrectly classified as '
                 'non-polyamines.\n'
                 '\n'
                 '2. **Amino Group Patterns**: The patterns used to identify '
                 'amino groups were too restrictive. For example, the pattern '
                 '`[NX3;H0]` only matches tertiary amines, but other forms of '
                 'tertiary amines (e.g., `[N+](=O)[O-]`) were not considered.\n'
                 '\n'
                 '3. **Complex Structures**: Some complex structures with '
                 'multiple amino groups were incorrectly classified, possibly '
                 'due to the way the patterns were matched or the way the '
                 'SMILES string was parsed.\n'
                 '\n'
                 '### Improvements:\n'
                 '- **Include Charged Amino Groups**: Add patterns to '
                 'recognize charged amino groups like `[NH3+]` and `[NH2+]`.\n'
                 '- **Expand Amino Group Patterns**: Include more variations '
                 'of amino groups, such as those in nitro groups or other '
                 'nitrogen-containing functional groups.\n'
                 '- **Handle Complex Structures**: Ensure that the code can '
                 'handle complex structures by using more robust pattern '
                 'matching.\n'
                 '\n'
                 'Here’s the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 131,
    'num_false_positives': 100,
    'num_true_negatives': 138,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.5670995670995671,
    'recall': 0.9849624060150376,
    'f1': 0.7197802197802199,
    'accuracy': 0.7250673854447439,
    'negative_predictive_value': 0.9857142857142858}