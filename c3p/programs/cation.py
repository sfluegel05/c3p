"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more 
    elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate formal charges on each atom
    positive_charges = 0
    positive_atoms = []
    
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge > 0:
            positive_charges += formal_charge
            positive_atoms.append(f"{atom.GetSymbol()}{atom.GetIdx()+1}")
    
    # Look for specific cationic patterns
    cationic_patterns = [
        Chem.MolFromSmarts('[N+]'), # Any positively charged nitrogen
        Chem.MolFromSmarts('[n+]'), # Any positively charged aromatic nitrogen
        Chem.MolFromSmarts('[P+]'), # Phosphonium
        Chem.MolFromSmarts('[S+]'), # Sulfonium
        Chem.MolFromSmarts('[O+]'), # Oxonium
        Chem.MolFromSmarts('[#7+]'), # Any positively charged nitrogen (alternative notation)
        Chem.MolFromSmarts('[NH4+]'), # Ammonium
        Chem.MolFromSmarts('[NH3+]'), # Primary ammonium
        Chem.MolFromSmarts('[NH2+]'), # Secondary ammonium
        Chem.MolFromSmarts('[NH+]'),  # Tertiary ammonium
        Chem.MolFromSmarts('[N+](C)(C)(C)'), # Quaternary ammonium
        Chem.MolFromSmarts('[Li+,Na+,K+,Rb+,Cs+,Fr+]'), # Alkali metals
        Chem.MolFromSmarts('[Be+2,Mg+2,Ca+2,Sr+2,Ba+2,Ra+2]'), # Alkaline earth metals
        Chem.MolFromSmarts('[#+1,#+2,#+3,#+4]') # Any atom with positive charge
    ]
    
    # If no explicit positive charges found, check for matches to cationic patterns
    if positive_charges == 0:
        for pattern in cationic_patterns:
            if pattern is not None and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    atom_idx = match[0]
                    atom = mol.GetAtomWithIdx(atom_idx)
                    positive_atoms.append(f"{atom.GetSymbol()}{atom_idx+1}")
                positive_charges += len(matches)
                
    # If still no positive charges found
    if positive_charges == 0:
        return False, "No positive charges or cationic groups found"
    
    # Success case - report all positive charges found
    charge_locations = ", ".join(positive_atoms)
    charge_str = f"{positive_charges}+" if positive_charges > 1 else "1+"
    return True, f"Found {charge_str} charge with cationic centers on: {charge_locations}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36916',
                          'name': 'cation',
                          'definition': 'A monoatomic or polyatomic species '
                                        'having one or more elementary charges '
                                        'of the proton.',
                          'parents': ['CHEBI:24870'],
                          'xrefs': ['KEGG:C01373'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7572815533980582 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[NH3+]C[C@H]1O[C@H](O[C@@H]2[C@@H]([NH3+])C[C@@H]([NH3+])[C@H](O)[C@H]2O[C@@H]2O[C@H](COP([O-])([O-])=O)[C@@H](O)[C@H]2O)[C@H]([NH3+])[C@@H](O)[C@@H]1O '
               "NAME: 5''-phosphoribostamycin(2+) REASON: CORRECT Found net "
               'positive charge (2+) with positive charges on: N1, N9, N12, '
               'N31\n'
               ' * SMILES: CCCCCCCCCC[C@@H](O)[C@@H](O)[C@@H]([NH3+])CO NAME: '
               'tetradecaphytosphingosine(1+) REASON: CORRECT Found net '
               'positive charge (1+) with positive charges on: N16\n'
               ' * SMILES: '
               '[H]C(=C([H])c1oc2ccccc2[n+]1CCCCCCCCCCCCCCCCCC)C([H])=C1Oc2ccccc2N1CCCCCCCCCCCCCCCCCC '
               'NAME: diOC18(3)(1+) REASON: CORRECT Found net positive charge '
               '(1+) with positive charges on: N11\n'
               ' * SMILES: '
               'CC[C@@]1(O)C[C@H](O[C@H]2C[C@@H]([C@H](O[C@H]3C[C@H](O)[C@H](O[C@H]4CC[C@H](O)[C@H](C)O4)[C@H](C)O3)[C@H](C)O2)[NH+](C)C)c2c(O)c3C(=O)c4c(O)cccc4C(=O)c3cc2[C@H]1C(=O)OC '
               'NAME: aclacinomycin N(1+) REASON: CORRECT Found net positive '
               'charge (1+) with positive charges on: N33\n'
               ' * SMILES: [H][C@@]12Cc3ccccc3C[N+]1(C)CCc1ccccc21 NAME: '
               'N-methyl-alpha-berbine REASON: CORRECT Found net positive '
               'charge (1+) with positive charges on: N10\n'
               ' * SMILES: '
               'C1(O)=C2C(C(C=3C(=CC=CC3O)C2=O)=O)=C(O)C4=C1C[C@](C[C@@H]4O[C@H]5C[C@@H]([C@@H]([C@@H](O5)C)O)[NH3+])([C@@H](C)O)O '
               'NAME: (13R)-13-dihydrocarminomycin(1+) REASON: CORRECT Found '
               'net positive charge (1+) with positive charges on: N33\n'
               ' * SMILES: '
               'NC(=[NH2+])N[C@H]1[C@H](O)[C@@H](O)[C@H](O)C(=O)[C@@H]1O NAME: '
               '1D-1-guanidiniumyl-1-deoxy-3-dehydro-scyllo-inositol(1+) '
               'REASON: CORRECT Found net positive charge (1+) with positive '
               'charges on: N3\n'
               ' * SMILES: '
               'CSCC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@@H]([NH3+])CC1=CC=CC=C1)C(N)=O '
               'NAME: L-Phe-Gly-L-Leu-L-Met-NH2(1+) REASON: CORRECT Found net '
               'positive charge (1+) with positive charges on: N22\n'
               ' * SMILES: CC1(C)CC(=O)N(CCCC[NH+]2CCN(CC2)C2=NC=CC=N2)C(=O)C1 '
               'NAME: gepirone(1+) REASON: CORRECT Found net positive charge '
               '(1+) with positive charges on: N12\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC=2C3=CC(O[C@@H]4O[C@H]([C@H](O)[C@@H](O)[C@H]4O)C)=C([O+]=C3C=C(O)C2)C=5C=C(OC)C(O)=C(O)C5)CO '
               'NAME: Petunidin 3-rhamnoside 5-glucoside REASON: CORRECT Found '
               'net positive charge (1+) with positive charges on: O27\n'
               ' * SMILES: CCC(C)CCCCCCCC\\C=C\\[C@@H](O)[C@@H]([NH3+])CO '
               'NAME: 14-methylhexadecasphingosine(1+) REASON: CORRECT Found '
               'net positive charge (1+) with positive charges on: N18\n'
               ' * SMILES: '
               'CC(=O)N[C@H]1C(O)O[C@H](COP(O)(=O)OCC[N+](C)(C)C)[C@H](O)[C@@H]1O '
               'NAME: N-acetyl-D-galactosamine-6-phosphocholine REASON: '
               'CORRECT Found net positive charge (1+) with positive charges '
               'on: N18\n'
               ' * SMILES: '
               'C1(=CC=C2C(=C1)C(=C3C(=N2)C=4N(C3)C(C5=C(C4)[C@](C(OC5)=O)(CC)O)=O)CC)OC(=O)N6CCC(CC6)[NH+]7CCCCC7 '
               'NAME: irinotecan(1+) REASON: CORRECT Found net positive charge '
               '(1+) with positive charges on: N38\n'
               ' * SMILES: '
               '[H][C@]12[C@@H](O)C[C@H]([NH3+])C(=O)N1CC1=C(C=CC(O)=C1O)[C@H]2O '
               'NAME: fumisoquin A(1+) REASON: CORRECT Found net positive '
               'charge (1+) with positive charges on: N6\n'
               ' * SMILES: '
               'COc1cc(C[C@H]2c3c(CC[N@+]2(C)CCCOC(=O)CCC(=O)OCCC[N@+]2(C)CCc4cc(OC)c(OC)c(OC)c4[C@H]2Cc2cc(OC)c(OC)c(OC)c2)cc(OC)c(OC)c3OC)cc(OC)c1OC '
               'NAME: meso-doxacurium REASON: CORRECT Found net positive '
               'charge (2+) with positive charges on: N12, N28\n'
               ' * SMILES: '
               '[C@@H]1(N(C([C@@H](NC([C@@H](NC([C@@H](NC([C@H](CC2=CC=C(O)C=C2)NC([C@@H](NC([C@H]3NC(=O)CC3)=O)CC(C)C)=O)=O)CCC([O-])=O)=O)CC(N)=O)=O)CCCC[NH3+])=O)CCC1)C(N[C@H](C(N[C@H](C(N4[C@H](C(N[C@@H](CC=5C=CC(=CC5)O)C(=O)[O-])=O)CCC4)=O)CCCNC(N)=[NH2+])=O)CCCNC(N)=[NH2+])=O '
               'NAME: neurotensin(1-11)(1+) REASON: CORRECT Found net positive '
               'charge (1+) with positive charges on: N55, N94, N102\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)\\C=C\\c2ccc(O)c(O)c2)O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)cc2)[C@H](O)[C@H]1O '
               'NAME: Pelargonidin 3-O-(6-caffeoyl-beta-D-glucoside) REASON: '
               'CORRECT Found net positive charge (1+) with positive charges '
               'on: O31\n'
               ' * SMILES: '
               'CCOC(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC(Cl)=C(OCC[NH+]2CC[NH+](C)CC2)C(Cl)=C1O '
               'NAME: JTE-607(2+) REASON: CORRECT Found net positive charge '
               '(2+) with positive charges on: N25, N28\n'
               ' * SMILES: C[N+]1(CCCCC[N+]2(C)CCCC2)CCCC1 NAME: pentolinium '
               'ion REASON: CORRECT Found net positive charge (2+) with '
               'positive charges on: N2, N8\n'
               ' * SMILES: [H][N+]([H])([H])[H] NAME: ammonium REASON: CORRECT '
               'Found net positive charge (1+) with positive charges on: N1\n'
               ' * SMILES: '
               '[O-][N+](=O)c1ccc(cc1)-n1nc(n[n+]1-c1ccc(I)cc1)-c1ccccc1 NAME: '
               'iodonitrotetrazolium REASON: CORRECT Found net positive charge '
               '(1+) with positive charges on: N2, N14\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)CC(O)=O)O[C@@H](Oc2cc(O)cc3[o+]c(c(O[C@@H]4O[C@H](COC(=O)\\C=C\\c5ccc(O)cc5)[C@@H](O)[C@H](O)[C@H]4O)cc23)-c2cc(O)c(O)c(O)c2)[C@H](O)[C@H]1O '
               'NAME: Malonylawobanin REASON: CORRECT Found net positive '
               'charge (1+) with positive charges on: O21\n'
               ' * SMILES: Cc1ncc(C[NH3+])c(N)n1 NAME: '
               '4-amino-5-ammoniomethyl-2-methylpyrimidine REASON: CORRECT '
               'Found net positive charge (1+) with positive charges on: N7\n'
               ' * SMILES: [La++] NAME: lanthanum(2+) REASON: CORRECT Found '
               'net positive charge (2+) with positive charges on: La1\n'
               ' * SMILES: C[n+]1ccc2c(c1)[nH]c1ccc(O)cc21 NAME: '
               '2-methyl-beta-carbolinium-6-ol REASON: CORRECT Found net '
               'positive charge (1+) with positive charges on: N2\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C([C@@](COC(CCCCCCC/C=C\\CCCCCC)=O)(OC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: PC(16:1(9Z)/20:3(5Z,8Z,11Z)) REASON: MISSED Net charge '
               '(0) is not positive\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCC)([O-])=O '
               'NAME: PC(15:0/18:4(6Z,9Z,12Z,15Z)) REASON: MISSED Net charge '
               '(0) is not positive\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)([O-])=O '
               'NAME: PC(22:1(11Z)/0:0) REASON: MISSED Net charge (0) is not '
               'positive\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCC[C@H]([C@H](COP(=O)([O-])OCC[N+](C)(C)C)NC(CCCCCCCCCCCCCCCCC)=O)O '
               'NAME: SM(d18:0/18:0) REASON: MISSED Net charge (0) is not '
               'positive\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCC)([O-])=O '
               'NAME: PC(14:0/20:5(5Z,8Z,11Z,14Z,17Z)) REASON: MISSED Net '
               'charge (0) is not positive\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)C)COC(=O)CCCCCCCCCCCCC)([O-])=O '
               'NAME: PC(14:0/2:0) REASON: MISSED Net charge (0) is not '
               'positive\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C/C=C\\CCCCCC)COC(=O)CCCCCCCCCCCCCCC)([O-])=O '
               'NAME: PC(16:0/18:2(9E,11Z)) REASON: MISSED Net charge (0) is '
               'not positive\n'
               ' * SMILES: '
               'C([C@@](COC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: PC(18:1(9Z)/20:2(11Z,14Z)) REASON: MISSED Net charge (0) '
               'is not positive\n'
               ' * SMILES: '
               'C([C@@](COC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCC/C=C\\CCCCCC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: PC(22:5(4Z,7Z,10Z,13Z,16Z)/16:1(9Z)) REASON: MISSED Net '
               'charge (0) is not positive\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCCCC '
               'NAME: 1,2-distearoyl-sn-glycero-3-phosphocholine REASON: '
               'MISSED Net charge (0) is not positive\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No positive charges or '
                                               'cationic groups found'}],
    'sample_false_negatives': [   {   'smiles': 'Cc1ccc(cc1)S(O)(=O)=O.N[C@@H]1CCN(C1)c1nc2n(cc(C(O)=O)c(=O)c2cc1F)-c1ccc(F)cc1F',
                                      'name': '(R)-tosufloxacin tosylate',
                                      'reason': 'No positive charges or '
                                                'cationic groups found'},
                                  {   'smiles': 'CO[C@H]1\\C=C\\O[C@@]2(C)OC3=C(C2=O)C2=C(C(O)=C3C)C([O-])=C(NC(=O)\\C(C)=C/C=C/[C@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@@H](C)[C@H](OC(C)=O)[C@@H]1C)\\C(=C/N(O)N1CCN(C)CC1)C2=O',
                                      'name': "2'-hydroxyrifampicin(1-)",
                                      'reason': 'No positive charges or '
                                                'cationic groups found'},
                                  {   'smiles': 'C1(O)=C(C=C(C(=C1)NC=2N=CN=C3C2C=CC(=C3)OCC4=CC=CC=C4)F)C',
                                      'name': 'ZM 323881',
                                      'reason': 'No positive charges or '
                                                'cationic groups found'},
                                  {   'smiles': '[H][C+]([H])([H])[H]',
                                      'name': 'methaniumyl',
                                      'reason': 'Invalid SMILES string'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                      'name': 'PG(20:0/21:0)',
                                      'reason': 'No positive charges or '
                                                'cationic groups found'},
                                  {   'smiles': 'COC1=C(Br)[C@H](O)[C@]2(CC(=NO2)C(=O)NCCCCNC(=O)C2=NO[C@]3(C2)C=C(Br)C(=O)[C@H](Br)[C@H]3O)C=C1Br',
                                      'name': '(+)-aplysinillin',
                                      'reason': 'No positive charges or '
                                                'cationic groups found'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O',
                                      'name': 'PA(20:2(11Z,14Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                      'reason': 'No positive charges or '
                                                'cationic groups found'}],
    'attempt': 1,
    'reasoning': 'Looking at the outcomes, I notice that the main issue is '
                 'that the program is failing to identify some phosphocholine '
                 '(PC) lipids and similar molecules as cations. These '
                 'molecules typically have a positively charged quaternary '
                 'ammonium group [N+(CH3)3] and a negatively charged phosphate '
                 'group [PO4-], giving them a net neutral charge overall.\n'
                 '\n'
                 'The key insight is that the previous code was too strict in '
                 'requiring a net positive charge. According to the '
                 'definition, a cation is "a monoatomic or polyatomic species '
                 'having one or more elementary charges of the proton." This '
                 'means that a molecule should be classified as a cation if it '
                 'contains ANY positively charged atoms, regardless of whether '
                 'there are also negative charges that make the overall '
                 'molecule neutral.\n'
                 '\n'
                 "Here's how we can improve the code:\n"
                 '\n'
                 '1. Remove the net charge requirement - we only need to check '
                 'for presence of positive charges\n'
                 '2. Keep track of the number of positive charges for '
                 'reporting\n'
                 '3. Add specific pattern matching for common cationic groups '
                 'like [N+] and quaternary ammonium\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2705,
    'num_false_positives': 100,
    'num_true_negatives': 2341,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.964349376114082,
    'recall': 0.9974188790560472,
    'f1': 0.980605401486315,
    'accuracy': 0.9792353968562003,
    'negative_predictive_value': 0.9970187393526405}