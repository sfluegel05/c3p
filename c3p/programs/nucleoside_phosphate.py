"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate has a nucleobase, a sugar moiety, and one or more phosphate groups attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extended nucleobase patterns (covering purine and pyrimidine with common modifications)
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1[nH]cnc2c1ncnc2"),   # Purines (adenine-like)
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),      # Purines (guaninie-like)
        Chem.MolFromSmarts("c1ccncn1"),           # Pyrimidine-like
        Chem.MolFromSmarts("c1nc[nH]c1=O"),       # Pyrimidone-like
        Chem.MolFromSmarts("n1c(=O)c[nH]cn1"),    # Cytosine-like
        Chem.MolFromSmarts("c1ccc(nc1)=O")        # Modified pyrimidines
    ]
    
    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found"
    
    # Refined sugar moiety detection
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(CO)C(O)C1"),  # Ribose / Deoxyribose and basic modifications considered
        Chem.MolFromSmarts("C1OCC(O)C(O)C1"),  # Alternate sugar backbone for modifications
        Chem.MolFromSmarts("C1OCC(C)C1")       # Tetrahydrofuran-like (ribose precursor)
    ]
    
    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_found:
        return False, "No sugar moiety found"
    
    # Flexible phosphate group detection, allowing for mono- and polyphosphate
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)(O)O"),  # Phosphate
        Chem.MolFromSmarts("P(=O)(O)(O)O"),  # Diphosphate link possibility, if misunderstood as poly
        Chem.MolFromSmarts("OP(O)(=O)OP(=O)(O)O") # Triphosphate (ATP/GTP like)
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns):
        return False, "No phosphate group found"

    return True, "Contains nucleobase, sugar moiety, and phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25608',
                          'name': 'nucleoside phosphate',
                          'definition': 'A nucleobase-containing molecular '
                                        'entity that is a nucleoside in which '
                                        'one or more of the sugar hydroxy '
                                        'groups has been converted into a '
                                        'mono- or poly-phosphate. The term '
                                        'includes both nucleotides and '
                                        'non-nucleotide nucleoside phosphates.',
                          'parents': [   'CHEBI:25703',
                                         'CHEBI:37734',
                                         'CHEBI:61120'],
                          'xrefs': ['KEGG:C01329'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.8571428571428571 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CC(C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)=C1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C '
               'NAME: 20-methyl-3-oxopregna-4,17-dien-21-oyl-CoA REASON: '
               'CORRECT Contains nucleobase, sugar moiety, and phosphate '
               'group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,13Z,16Z,19Z,22Z)-3-hydroxyoctacosatetraenoyl-CoA '
               'REASON: CORRECT Contains nucleobase, sugar moiety, and '
               'phosphate group\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCc1ccc(O)c(O)c1 '
               'NAME: dihydrocaffeoyl-CoA REASON: CORRECT Contains nucleobase, '
               'sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,11Z)-3-hydroxyoctadecenoyl-CoA REASON: CORRECT '
               'Contains nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)c1ccc(F)cc1 '
               'NAME: 4-fluorobenzoyl-CoA REASON: CORRECT Contains nucleobase, '
               'sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'C[C@@H](O)[C@H](NC(=O)Nc1nc(S)nc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O)C(O)=O '
               "NAME: 2-thio-N(6)-L-threonylcarbamoyladenine 5'-monophosphate "
               'REASON: CORRECT Contains nucleobase, sugar moiety, and '
               'phosphate group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (4Z,7Z,10Z,13Z,16Z)-docosapentaenoyl-CoA REASON: CORRECT '
               'Contains nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COP(O)(=O)OP(O)(O)=O)O[C@H]([C@@H]1O)n1cnc2c1[nH]c(=O)[nH]c2=O '
               'NAME: XDP REASON: CORRECT Contains nucleobase, sugar moiety, '
               'and phosphate group\n'
               ' * SMILES: '
               '[C@@H]1([C@@H](C[C@H]([C@@H](O1)C)O)O)OCCCCCCCCCCCCCCCCCCC(=O)SCCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O '
               'NAME: oscr#34-CoA REASON: CORRECT Contains nucleobase, sugar '
               'moiety, and phosphate group\n'
               ' * SMILES: '
               'CCC(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12 '
               'NAME: propanoyl-AMP REASON: CORRECT Contains nucleobase, sugar '
               'moiety, and phosphate group\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,13Z)-3-hydroxyicosenoyl-CoA REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/CCCC\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (5Z,11Z,14Z,17Z)-icosatetraenoyl-CoA REASON: CORRECT '
               'Contains nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (13Z,16Z,19Z,22Z,25Z)-octacosapentaenoyl-CoA REASON: '
               'CORRECT Contains nucleobase, sugar moiety, and phosphate '
               'group\n'
               ' * SMILES: '
               'S(C(CCCCCCC/C=C\\C[C@@H](CCCCCC)O)=O)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]1O[C@@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H]([C@@H]1OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O '
               'NAME: (9Z,12R)-12-hydroxyoctadec-9-enoyl-CoA REASON: CORRECT '
               'Contains nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'S(CCNC(=O)CCNC(=O)C(O)C(COP(OP(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)[C@H](O)[C@@H]1OP(O)(O)=O)(O)=O)(O)=O)(C)C)C(=O)C=4C(=CC=CC4)C '
               'NAME: Salicyl-CoA REASON: CORRECT Contains nucleobase, sugar '
               'moiety, and phosphate group\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCSC(=O)C(=S)CCC(O)=O '
               'NAME: 4-carboxy-2-thioxobutanoyl-CoA REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (8Z,11Z,14Z,17Z)-icosatetraenoyl-CoA REASON: CORRECT '
               'Contains nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC(=O)CCCCCCCCC(O)=O '
               'NAME: 3-oxododecanedioyl-CoA REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCC/C=C\\CCCCCCCC)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (7Z)-hexadecenoyl-CoA REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)\\C=C/CCC(O)=O '
               'NAME: (Z)-2,3-dehydroadipoyl-CoA REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N '
               'NAME: oleoyl-CoA REASON: CORRECT Contains nucleobase, sugar '
               'moiety, and phosphate group\n'
               ' * SMILES: '
               'Nc1ncnc2n([C@@H]3O[C@@H]4COP(O)(=O)O[C@H]4[C@H]3O)c(Sc3ccc(Cl)cc3)nc12 '
               'NAME: 8-(4-chlorophenylthio)-cAMP REASON: CORRECT Contains '
               'nucleobase, sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'N(C1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(OP(=O)(OP(=O)([O-])[O-])[O-])[O-])[C@@H](O)[C@H]3O)C '
               'NAME: N(6)-methyl-ATP(4-) REASON: CORRECT Contains nucleobase, '
               'sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(O)(=O)OC(=O)CCCC[C@@H]2CCSS2)[C@@H](O)[C@H]1O '
               'NAME: (R)-lipoyl-GMP REASON: CORRECT Contains nucleobase, '
               'sugar moiety, and phosphate group\n'
               ' * SMILES: '
               'Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(O)(=O)OP(O)(O)=O)O1 '
               'NAME: dADP REASON: CORRECT Contains nucleobase, sugar moiety, '
               'and phosphate group\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'Cc1cc2[N]c3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)c2cc1C '
               'NAME: FMNH(.) REASON: MISSED No nucleobase found\n'
               ' * SMILES: '
               'C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O '
               "NAME: 8-formyl-8-demethylriboflavin 5'-phosphate REASON: "
               'MISSED No nucleobase found\n'
               ' * SMILES: '
               'Nc1ccn([C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H]2O)c(=O)n1 '
               "NAME: 3'-CMP REASON: MISSED No nucleobase found\n"
               ' * SMILES: '
               'C[C@H]1O[C@@H](C[C@@](C)(O)C1=O)OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C[C@@H]1O)n1cc(C)c(=O)[nH]c1=O '
               'NAME: dTDP-3-methyl-4-oxo-2,6-dideoxy-D-glucose REASON: MISSED '
               'No nucleobase found\n'
               ' * SMILES: '
               'N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1cnc2c(N)ncnc12 '
               "NAME: N(6),N(6),O-tridemethylpuromycin 5'-phosphate REASON: "
               'MISSED No sugar moiety found\n'
               ' * SMILES: '
               'C[C@H]1O[C@@H](C[C@@H](N)[C@@H]1O)OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C[C@@H]1O)n1cc(C)c(=O)[nH]c1=O '
               'NAME: dTDP-2,6-dideoxy-D-kanosamine REASON: MISSED No '
               'nucleobase found\n'
               ' * SMILES: '
               'C[C@@H]1CC(=O)[C@@H](O)[C@@H](OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cc(C)c(=O)[nH]c2=O)O1 '
               'NAME: TDP-actinospectose REASON: MISSED No nucleobase found\n'
               ' * SMILES: '
               'C1=CC(=NC(N1[C@@H]2O[C@H](COP(O[C@]3(C[C@H](O)[C@@H](O)[C@@](O3)([C@@H]([C@@H](CO)O)O)[H])C(=O)O)(=O)O)[C@H]([C@H]2O)O)=O)N '
               'NAME: CMP-3-deoxy-D-glycero-beta-D-galacto-nonulosonic acid '
               'REASON: MISSED No nucleobase found\n'
               ' * SMILES: '
               'NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O '
               'NAME: NMNH REASON: MISSED No nucleobase found\n'
               ' * SMILES: '
               'Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)c2cc1C '
               'NAME: FMNH2 REASON: MISSED No nucleobase found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No nucleobase found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No nucleobase found'}],
    'sample_false_negatives': [   {   'smiles': 'Cc1cc2[N]c3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)c2cc1C',
                                      'name': 'FMNH(.)',
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'C(N1C=2C(=NC3=C1C=C(C(=C3)C)C(=O)[H])C(NC(N2)=O)=O)[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)O',
                                      'name': '8-formyl-8-demethylriboflavin '
                                              "5'-phosphate",
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1cnc2c(N)ncnc12',
                                      'name': 'N(6),N(6),O-tridemethylpuromycin '
                                              "5'-phosphate",
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O',
                                      'name': 'NMNH',
                                      'reason': 'No nucleobase found'},
                                  {   'smiles': 'Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)c2cc1C',
                                      'name': 'FMNH2',
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'N1=C(C)C(NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)=O',
                                      'name': "6-azathymidine 5'-monophosphate",
                                      'reason': 'No nucleobase found'},
                                  {   'smiles': 'P(=O)(OC[C@@H](O)[C@@H](O)[C@@H](O)CN1C=2NC(=O)NC(C2[N+]3=CCC(C=4C3=C1C=C(C)C4C)(C)C)=O)(O)O',
                                      'name': 'Prenylated flavin '
                                              'mononucleotide',
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'O(C[C@H]1O[C@@H]([N+]2=CC3=C4C(=C2)C(=[S+][Ni-]4SC3=O)O)[C@@H]([C@@H]1O)O)P(O)(=O)O',
                                      'name': 'Ni(II)-pyridinium-3,5-bisthiocarboxylic '
                                              'acid mononucleotide',
                                      'reason': 'No nucleobase found'},
                                  {   'smiles': 'P(OC[C@H]1O[C@@H](N2C3=NC=NC(N)=C3N=C2)CC1)(OP(OP(O)(O)=O)(O)=O)(O)=O',
                                      'name': "2',3'-Dideoxyadenosine-5-triphosphate",
                                      'reason': 'No sugar moiety found'},
                                  {   'smiles': 'Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C[C@H]1O',
                                      'name': 'cordycepin triphosphate',
                                      'reason': 'No sugar moiety found'}],
    'attempt': 2,
    'reasoning': 'To address the challenges observed in the previous attempts, '
                 'we will refine the detection logic for identifying '
                 'nucleoside phosphates. The prior implementation seemed to '
                 'have issues recognizing certain nucleobase and sugar moiety '
                 'patterns correctly. The following strategy will be used for '
                 'improvement:\n'
                 '\n'
                 '1. **Nucleobase Recognition**: The nucleobase patterns need '
                 'to be expanded to cover common nucleobases and their '
                 'variants. Inclusion of modified bases would be helpful in '
                 'handling examples like FMNH and other modified '
                 'environments.\n'
                 '\n'
                 '2. **Sugar Moiety Detection**: We must improve the detection '
                 'of sugar moieties by recognizing both ribose and '
                 'deoxyribose, with considerations for common modifications.\n'
                 '\n'
                 '3. **Phosphate Group Detection**: Ensure that the phosphate '
                 'group detection is flexible enough to detect both mono- and '
                 'poly-phosphate groups. \n'
                 '\n'
                 '4. **Connectivity Check**: Beyond detection of '
                 'substructures, a check to ensure that the nucleobase, sugar, '
                 'and phosphate groups are connected correctly to form a '
                 'nucleoside phosphate.\n'
                 '\n'
                 "Here's the refined program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1264,
    'num_false_positives': 100,
    'num_true_negatives': 8709,
    'num_false_negatives': 38,
    'num_negatives': None,
    'precision': 0.9266862170087976,
    'recall': 0.9708141321044547,
    'f1': 0.9482370592648162,
    'accuracy': 0.986351498368114,
    'negative_predictive_value': 0.9956556533668687}