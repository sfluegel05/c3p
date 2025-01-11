"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an amine and ester group attached to the same carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern 1: Basic alpha-amino acid ester pattern
    # Matches both primary and secondary amines
    pattern1 = Chem.MolFromSmarts("[NX3;!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[C,H]")
    
    # Pattern 2: N-substituted amino acid ester pattern (including tertiary amines)
    pattern2 = Chem.MolFromSmarts("[#7;!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[C,H]")
    
    # Pattern 3: Cyclic amino acid ester pattern (like proline derivatives)
    pattern3 = Chem.MolFromSmarts("[#7;R][CH1,CH2;R][C](=O)O[C,H]")
    
    # Pattern 4: Complex ring systems with embedded amino acid ester
    pattern4 = Chem.MolFromSmarts("[#7;R][CH1,CH2][C](=O)O[C,H]")
    
    patterns = [
        (pattern1, "basic alpha-amino acid ester"),
        (pattern2, "N-substituted amino acid ester"),
        (pattern3, "cyclic amino acid ester"),
        (pattern4, "ring-containing amino acid ester")
    ]
    
    for pattern, pattern_name in patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                n_atom = mol.GetAtomWithIdx(match[0])
                c_atom = mol.GetAtomWithIdx(match[1])
                
                # Verify nitrogen properties
                if n_atom.GetAtomicNum() != 7:  # Must be nitrogen
                    continue
                    
                # Check if carbon is really alpha to both N and C(=O)
                has_n = False
                has_ester = False
                for neighbor in c_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen
                        has_n = True
                    elif neighbor.GetAtomicNum() == 6:  # Carbon
                        for bond in neighbor.GetBonds():
                            other_atom = bond.GetOtherAtom(neighbor)
                            if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_ester = True
                
                if has_n and has_ester:
                    # Additional check for valid ester group
                    ester_env = Chem.MolFromSmarts("[C](=O)O[C,H]")
                    if mol.HasSubstructMatch(ester_env):
                        return True, f"Contains {pattern_name} pattern"
    
    return False, "No valid alpha-amino acid ester pattern found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46874',
                          'name': 'alpha-amino acid ester',
                          'definition': 'The amino acid ester derivative '
                                        'obtained the formal condensation of '
                                        'an alpha-amino acid with an alcohol.',
                          'parents': ['CHEBI:46668'],
                          'xrefs': ['KEGG:C03317'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5614035087719299 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'COC(=O)[C@@H]1[C@H]([C@@H]2CN3C(=O)C=CC=C3[C@@H]2N1CC4CCOCC4)CO '
               'NAME: methyl '
               '(2S,3S,3aR,9bR)-3-(hydroxymethyl)-1-(oxan-4-ylmethyl)-6-oxo-3,3a,4,9b-tetrahydro-2H-pyrrolo[2,3-a]indolizine-2-carboxylate '
               'REASON: CORRECT Contains alpha-amino acid ester pattern\n'
               ' * SMILES: ClC=1C([C@H](N2CCC(S)=C(C2)CC(O)=O)C(OC)=O)=CC=CC1 '
               'NAME: '
               '{1-[(1s)-1-(2-chlorophenyl)-2-methoxy-2-oxoethyl]-4-sulfanyl-1,2,5,6-tetrahydro-3-pyridinyl}acetic '
               'acid REASON: CORRECT Contains alpha-amino acid ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC([C@H](CCCCN)N)=O)[C@H]3O '
               "NAME: 3'-L-lysyl-AMP REASON: CORRECT Contains alpha-amino acid "
               'ester pattern\n'
               ' * SMILES: C=1(C[C@@H](C(OC)=O)N)C=CC(O)=CC1 NAME: methyl '
               'L-tyrosinate REASON: CORRECT Contains alpha-amino acid ester '
               'pattern\n'
               ' * SMILES: '
               'COC(=O)[C@H]1[C@@H]([C@H]2CN3C(=CC=C(C3=O)C4=CCCCC4)[C@H]2N1CC5=CC6=C(C=C5)OCO6)CO '
               'NAME: '
               '(2R,3R,3aS,9bS)-1-(1,3-benzodioxol-5-ylmethyl)-7-(1-cyclohexenyl)-3-(hydroxymethyl)-6-oxo-3,3a,4,9b-tetrahydro-2H-pyrrolo[2,3-a]indolizine-2-carboxylic '
               'acid methyl ester REASON: CORRECT Contains alpha-amino acid '
               'ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CC(=O)O)[C@H]3O '
               "NAME: 3'-L-aspartyl-AMP REASON: CORRECT Contains alpha-amino "
               'acid ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)CN)[C@H]3O '
               "NAME: 3'-glycyl-AMP REASON: CORRECT Contains alpha-amino acid "
               'ester pattern\n'
               ' * SMILES: '
               '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COC(=O)CN(CC)CC '
               'NAME: hydrocortamate REASON: CORRECT Contains alpha-amino acid '
               'ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC([C@H](CC=4N=CNC4)N)=O)[C@H]3O '
               "NAME: 3'-L-histidyl-AMP REASON: CORRECT Contains alpha-amino "
               'acid ester pattern\n'
               ' * SMILES: COC(=O)C1CCC(CO)N1 NAME: methyl '
               '5-(hydroxymethyl)pyrrolidine-2-carboxylate REASON: CORRECT '
               'Contains alpha-amino acid ester pattern\n'
               ' * SMILES: C(C(OCCN(CC)CC)=O)(C1=CC=CC=C1)(N2CCCCC2)[H] NAME: '
               'bietamiverine REASON: CORRECT Contains alpha-amino acid ester '
               'pattern\n'
               ' * SMILES: '
               'O=C(OC1[C@@H](OC[C@@H]([C@@H]1OC)C)/C(=C/C=C/C=C/C)/C)CN NAME: '
               'Lanomycin REASON: CORRECT Contains alpha-amino acid ester '
               'pattern\n'
               ' * SMILES: O(CCCCCCCCCCCC)C(=O)C(N(C)C)C NAME: DDAIP REASON: '
               'CORRECT Contains alpha-amino acid ester pattern\n'
               ' * SMILES: '
               '[C@@H]([C@H](O)C1=CC=C(C=C1)S(C)(=O)=O)(NC(C(Cl)Cl)=O)COC(CN)=O '
               'NAME: thiamphenicol glycinate REASON: CORRECT Contains '
               'alpha-amino acid ester pattern\n'
               ' * SMILES: [C@@H](CCCNC(N[N+](=O)[O-])=N)(N)C(OC)=O NAME: '
               'N(gamma)-nitro-L-arginine methyl ester REASON: CORRECT '
               'Contains alpha-amino acid ester pattern\n'
               ' * SMILES: '
               'COC(=O)[C@@H]1[C@H]([C@@H]2CN3C(=CC=C(C3=O)C4=CN=CC=C4)[C@@H]2N1CC5=CC(=CC=C5)F)CO '
               'NAME: methyl '
               '(2S,3S,3aR,9bR)-1-[(3-fluorophenyl)methyl]-3-(hydroxymethyl)-6-oxo-7-pyridin-3-yl-3,3a,4,9b-tetrahydro-2H-pyrrolo[2,3-a]indolizine-2-carboxylate '
               'REASON: CORRECT Contains alpha-amino acid ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CS)[C@H]3O '
               "NAME: 3'-L-cysteinyl-AMP REASON: CORRECT Contains alpha-amino "
               'acid ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CC(=O)N)[C@H]3O '
               "NAME: 3'-L-asparaginyl-AMP REASON: CORRECT Contains "
               'alpha-amino acid ester pattern\n'
               ' * SMILES: O(C(=O)CN1CCNCC1)CC NAME: Ethyl piperazinoacetate '
               'REASON: CORRECT Contains alpha-amino acid ester pattern\n'
               ' * SMILES: C=1(C[C@@H](C(OCC)=O)N)C=CC(O)=CC1 NAME: ethyl '
               'L-tyrosinate REASON: CORRECT Contains alpha-amino acid ester '
               'pattern\n'
               ' * SMILES: CCN(CC)CCNC(C1=CC=CC=C1)C(=O)OCCC(C)C.Cl.Cl NAME: '
               'Acamylophenine dihydrochloride REASON: CORRECT Contains '
               'alpha-amino acid ester pattern\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CC=4C=CC(=CC4)O)[C@H]3O '
               "NAME: 3'-L-tyrosyl-AMP REASON: CORRECT Contains alpha-amino "
               'acid ester pattern\n'
               ' * SMILES: '
               'CC=CC1=CC=C2[C@H]3[C@@H](CN2C1=O)[C@@H]([C@H](N3CCC4=CC=CC=C4)C(=O)OC)CO '
               'NAME: LSM-14628 REASON: CORRECT Contains alpha-amino acid '
               'ester pattern\n'
               ' * SMILES: O(C(=O)CNCC(OC)=O)C NAME: Methyl '
               '2-[(2-methoxy-2-oxoethyl)amino]acetate REASON: CORRECT '
               'Contains alpha-amino acid ester pattern\n'
               ' * SMILES: OC(C(N)C(OC)=O)C NAME: Dl-threonine methyl ester '
               'REASON: CORRECT Contains alpha-amino acid ester pattern\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O=C1N2[C@]3(NC=4C(O)=CC=CC4[C@@]3(O)C[C@H]2C(=O)O[C@H]1CC5=CC=CC=C5)C(C=C)(C)C '
               'NAME: Shornephine A REASON: MISSED No alpha-amino acid ester '
               'pattern found\n'
               ' * SMILES: O=C1OCC2=CC=C(N2[C@H]1[C@H](CC)C)C=O NAME: '
               '(4S)-4-(2-methylbutyl)-3-oxo-3,4-dihydro-1H-pyrrolo[2,1-c][1,4]oxazine-6-carbaldehyde '
               'REASON: MISSED No valid alpha-amino acid ester pattern found\n'
               ' * SMILES: '
               'O=C1N2[C@@H]3NC=4C=CC=CC4[C@]3(CC=C(C)C)C[C@@H]2C(=O)O[C@@H]1CC(C)C '
               'NAME: Mollenine A REASON: MISSED No alpha-amino acid ester '
               'pattern found\n'
               ' * SMILES: '
               'C1=CC=CC=C1C[C@@H](C(=O)OCC)NC(=O)C2=CC(=C3C(=C2O)C(O[C@@H](C3)C)=O)Cl '
               'NAME: ochratoxin C REASON: MISSED No alpha-amino acid ester '
               'pattern found\n'
               ' * SMILES: [C@H]1(CCC(N1C)=O)C(=O)OC NAME: methyl '
               '1-methyl-5-oxo-L-prolinate REASON: MISSED No alpha-amino acid '
               'ester pattern found\n'
               ' * SMILES: '
               'COC1=C(C=C(C=C1)C(=O)COC(=O)CN2CC(=O)OC3=CC=CC=C32)Br NAME: '
               '2-(2-oxo-3H-1,4-benzoxazin-4-yl)acetic acid '
               '[2-(3-bromo-4-methoxyphenyl)-2-oxoethyl] ester REASON: MISSED '
               'No valid alpha-amino acid ester pattern found\n'
               ' * SMILES: '
               'BrC[C@H]1C2=C3C(N[C@](C3=O)(C(=O)OC)C)=C(O)C=C2N(C1)C(=O)C=4NC=5C(OC)=C(OC)C(=CC5C4)OC '
               'NAME: Duocarmycin B2 REASON: MISSED No alpha-amino acid ester '
               'pattern found\n'
               ' * SMILES: '
               'C1=CC=C2CN(C(C)=O)[C@H](C(OC3=CC(=CN=C3)Cl)=O)CC2=C1 NAME: '
               '(5-Chloropyridin-3-yl) '
               '(3S)-2-acetyl-3,4-dihydro-1H-isoquinoline-3-carboxylate '
               'REASON: MISSED No ester group found\n'
               ' * SMILES: CCCCOC(=O)C(C)NC1=NNC(=S)NC1=O NAME: '
               '2-[(5-oxo-3-sulfanylidene-2H-1,2,4-triazin-6-yl)amino]propanoic '
               'acid butyl ester REASON: MISSED No valid alpha-amino acid '
               'ester pattern found\n'
               ' * SMILES: '
               'CCOC1=C(C=C(C=C1)C=C2C(=O)N(C(=S)N2CC(=O)OC)C3=CC=C(C=C3)OC)OCC '
               'NAME: '
               '2-[5-[(3,4-diethoxyphenyl)methylidene]-3-(4-methoxyphenyl)-4-oxo-2-sulfanylidene-1-imidazolidinyl]acetic '
               'acid methyl ester REASON: MISSED No valid alpha-amino acid '
               'ester pattern found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No valid alpha-amino acid '
                                               'ester pattern found'}],
    'sample_false_negatives': [   {   'smiles': 'BrC[C@H]1C2=C3C(N[C@](C3=O)(C(=O)OC)C)=C(O)C=C2N(C1)C(=O)C=4NC=5C(OC)=C(OC)C(=CC5C4)OC',
                                      'name': 'Duocarmycin B2',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'C1=CC=C2CN(C(C)=O)[C@H](C(OC3=CC(=CN=C3)Cl)=O)CC2=C1',
                                      'name': '(5-Chloropyridin-3-yl) '
                                              '(3S)-2-acetyl-3,4-dihydro-1H-isoquinoline-3-carboxylate',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'ClCC1C2=C3C(NC(C3=O)(C(=O)OC)C)=C(O)C=C2N(C1)C(=O)C=4NC=5C(OC)=C(OC)C(=CC5C4)OC',
                                      'name': 'Duocarmycin C2',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C1N(NC2=C(OC)C=C(OC)C=C2)[C@@](OC)(C(=O)OC)CC1',
                                      'name': 'Peniproline A',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C1O[C@@H]([C@@H](N)C(=O)OC2=C1C=CC=C2NC=O)C',
                                      'name': 'Streptomyceamide B',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C(OC)[C@@H](N=[N+]([O-])/C=C\\CCCCCCCCCCCC(C)C)C',
                                      'name': 'Azodyrecin C',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C1N(C2(C(=O)OC2)C(C1(C)C)(O)C(CC(O)C3OC3C(O)CNC(=O)C(C(O)/C(=C\\C/C=C/C=C/C=C/C=4OC=NC4)/C)(C)C)C)C',
                                      'name': 'Oxazolepoxidomycin A',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'N[C@@H](CCCNC(N)=N)C(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C1N=CN=C2N',
                                      'name': 'L-arginyl-AMP',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C(OC)[C@@H](N=[N+]([O-])/C=C\\CCCCCCCCCC(C)C)C',
                                      'name': 'Azodyrecin A',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'},
                                  {   'smiles': 'O=C(OC)[C@@H](N=[N+]([O-])/C=C/CCCCCCCCCC(CC)C)C',
                                      'name': '1′-trans-Azodyrecin B',
                                      'reason': 'No valid alpha-amino acid '
                                                'ester pattern found'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Pattern Recognition Issues:\n'
                 '- The previous code missed several valid alpha-amino acid '
                 'esters because the SMARTS patterns were too restrictive\n'
                 '- It failed to recognize cyclic amino acid esters (like '
                 'proline derivatives)\n'
                 '- The nitrogen pattern `[NX3;!$(NC=O)]` was too strict and '
                 'missed some valid cases\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 '- Mollenine A and Shornephine A contain complex ring systems '
                 'with embedded amino acid ester patterns\n'
                 '- Methyl 1-methyl-5-oxo-L-prolinate was missed due to not '
                 'accounting for cyclic amino acids\n'
                 '- Some compounds with N-substituted amino acids were missed\n'
                 '\n'
                 '3. Improvements Needed:\n'
                 '- Broaden the SMARTS patterns to catch more variations\n'
                 '- Add specific patterns for cyclic amino acids\n'
                 '- Better handle N-substituted cases\n'
                 '- Account for more complex ring systems\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 212,
    'num_false_positives': 100,
    'num_true_negatives': 7477,
    'num_false_negatives': 24,
    'num_negatives': None,
    'precision': 0.6794871794871795,
    'recall': 0.8983050847457628,
    'f1': 0.7737226277372262,
    'accuracy': 0.9841290157429925,
    'negative_predictive_value': 0.9968004266097854}