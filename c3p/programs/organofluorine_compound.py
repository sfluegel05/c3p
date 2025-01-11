"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:24651 organofluorine compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as a compound containing at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon and fluorine atoms
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    has_fluorine = any(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms())
    
    if not has_carbon and not has_fluorine:
        return False, "Contains neither carbon nor fluorine atoms"
    if not has_carbon:
        return False, "Contains fluorine but no carbon atoms"
    if not has_fluorine:
        return False, "Contains no fluorine atoms"

    # Look for any carbon-fluorine bonds (both aromatic and aliphatic)
    cf_bonds = 0
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check if bond is between C and F
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or
            (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6)):
            cf_bonds += 1

    if cf_bonds == 0:
        return False, "No carbon-fluorine bonds found"

    # Count total fluorine atoms
    f_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9])
    
    return True, f"Contains {f_atoms} fluorine atom{'s' if f_atoms>1 else ''} ({cf_bonds} C-F bond{'s' if cf_bonds>1 else ''})"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37143',
                          'name': 'organofluorine compound',
                          'definition': 'An organofluorine compound is a '
                                        'compound containing at least one '
                                        'carbon-fluorine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:24062'],
                          'xrefs': ['MetaCyc:Fluorides'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7175141242937852 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: C(O)(=O)CCCCCCC/C=C\\C/C=C\\CCCCCF '
               'NAME: 18-fluoro-9Z,12Z-octadecadienoic acid REASON: CORRECT '
               'Contains 1 fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               'CN1CCN(Cc2ccc(NC(=O)c3ccc(C)c(c3)C#Cc3cnc4cccnn34)cc2C(F)(F)F)CC1 '
               'NAME: ponatinib REASON: CORRECT Contains 3 fluorine atoms (3 '
               'C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: FC(F)(F)C(F)(F)N1C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F '
               'NAME: perfluoro-N-ethylpiperidine REASON: CORRECT Contains 15 '
               'fluorine atoms (15 C-F bonds, 9 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'C[C@H]1C[C@H]2[C@@H]3CC[C@](OC(C)=O)(C(C)=O)[C@@]3(C)C[C@H](O)[C@]2(F)[C@@]2(C)C=CC(=O)C=C12 '
               'NAME: fluorometholone acetate REASON: CORRECT Contains 1 '
               'fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(N=CC=C2)O[C@H]1CN(C)CC3=CC(=CC=C3)C(F)(F)F)[C@H](C)CO '
               'NAME: '
               '(2R,3S)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl-[[3-(trifluoromethyl)phenyl]methyl]amino]methyl]-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-6-one '
               'REASON: CORRECT Contains 3 fluorine atoms (3 C-F bonds, 3 CF2 '
               'groups, 1 CF3 group)\n'
               ' * SMILES: '
               'FC1(F)C(F)(F)C(F)(F)C2(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C2(F)C1(F)F '
               'NAME: perfluorodecalin REASON: CORRECT Contains 18 fluorine '
               'atoms (18 C-F bonds, 8 CF2 groups)\n'
               ' * SMILES: [O-][N+](=O)CC(=O)c1ccccc1OC(F)(F)F NAME: '
               '2-nitro-1-(2-trifluoromethoxyphenyl)ethanone REASON: CORRECT '
               'Contains 3 fluorine atoms (3 C-F bonds, 3 CF2 groups, 1 CF3 '
               'group)\n'
               ' * SMILES: FC(F)(F)c1cccc(c1)C1=CCN(CCc2ccc3ccccc3c2)CC1 NAME: '
               'xaliproden REASON: CORRECT Contains 3 fluorine atoms (3 C-F '
               'bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: NCC(F)=O NAME: 2-aminoacetyl fluoride REASON: '
               'CORRECT Contains 1 fluorine atom (1 C-F bond)\n'
               ' * SMILES: FC(F)(F)C1=NC=CC=C1C(=O)NC1CCC1C1=C(Cl)C=C(Cl)C=C1 '
               'NAME: '
               'N-[2-(2,4-dichlorophenyl)cyclobutyl]-2-(trifluoromethyl)pyridine-3-carboxamide '
               'REASON: CORRECT Contains 3 fluorine atoms (3 C-F bonds, 3 CF2 '
               'groups, 1 CF3 group)\n'
               ' * SMILES: '
               'CO\\C=C(\\C(=O)OC)C1=C(OC2=NC=CC(=C2)C(F)(F)F)C=CC=C1 NAME: '
               'flupyroxystrobin REASON: CORRECT Contains 3 fluorine atoms (3 '
               'C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'C1(=NC(=NC(=N1)C2=NC(=CC=C2)C(F)(F)F)NCC(O)(C)C)NC=3C=CN=C(C3)C(F)(F)F '
               'NAME: enasidenib REASON: CORRECT Contains 6 fluorine atoms (6 '
               'C-F bonds, 6 CF2 groups, 2 CF3 groups)\n'
               ' * SMILES: '
               'C=1C(=C(C(=CC1C(F)(F)F)Cl)N2N=C(C(=O)C)C(=C2N)S(=O)C)Cl NAME: '
               'acetoprole REASON: CORRECT Contains 3 fluorine atoms (3 C-F '
               'bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'CN1[C@@H]2CC[C@H](O[C@@H]2COC3=C(C1=O)C=C(C=C3)NC(=O)C4=CC=C(C=C4)C(F)(F)F)CC(=O)OC '
               'NAME: '
               '2-[(2S,4aR,12aS)-5-methyl-6-oxo-8-[[oxo-[4-(trifluoromethyl)phenyl]methyl]amino]-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]acetic '
               'acid methyl ester REASON: CORRECT Contains 3 fluorine atoms (3 '
               'C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'C1(=NOC(=C1)COC=2C=CC=C(C2)C(F)(F)F)C(=O)N3CCC(CC3)C4=CC=NC=C4 '
               'NAME: dafadine D REASON: CORRECT Contains 3 fluorine atoms (3 '
               'C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F '
               'NAME: 2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9-hexadecafluorononanoic '
               'acid REASON: CORRECT Contains 16 fluorine atoms (16 C-F bonds, '
               '8 CF2 groups)\n'
               ' * SMILES: [H]C(=O)[C@H]([18F])[C@@H](O)[C@H](O)[C@H](O)CO '
               'NAME: 2-deoxy-2-((18)F)fluoro-aldehydo-D-glucose REASON: '
               'CORRECT Contains 1 fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               '[N@+]1([C@]2(C[C@H](OC(C(C3=CC=CC=C3)(C4=CC=CC=C4)O)=O)C[C@@]1(CC2)[H])[H])(CCF)C '
               'NAME: flutropium REASON: CORRECT Contains 1 fluorine atom (1 '
               'C-F bond)\n'
               ' * SMILES: '
               'CON=CNC(=O)c1ccc(cc1C)C1=NO[C@@](C1)(c1cc(Cl)cc(Cl)c1)C(F)(F)F '
               'NAME: (S)-fluxametamide REASON: CORRECT Contains 3 fluorine '
               'atoms (3 C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: F[C@@]12C3(C(=C[C@H](O)CC3)CCC1C=4CCC(C4CC2O)(C)C)C '
               'NAME: '
               '(7R,9bR)-9b-fluoro-1,1,9a-trimethyl-2H,3H,3bH,4H,5H,7H,8H,9H,10H,11H-cyclopenta[a]phenanthrene-7,10-diol '
               'REASON: CORRECT Contains 1 fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)CC2=C(C=CC(=C2)NC(=O)CC3=CC=CC=C3)O[C@H]1CN(C)CC4=CC=C(C=C4)C(F)(F)F)[C@@H](C)CO '
               'NAME: '
               'N-[(2R,3S)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl-[[4-(trifluoromethyl)phenyl]methyl]amino]methyl]-6-oxo-2,3,4,7-tetrahydro-1,5-benzoxazonin-9-yl]-2-phenylacetamide '
               'REASON: CORRECT Contains 3 fluorine atoms (3 C-F bonds, 3 CF2 '
               'groups, 1 CF3 group)\n'
               ' * SMILES: '
               'CCN1OC[C@@H](NC(=O)C2=C(C)C=C(C=C2)C2=NO[C@](C2)(C2=CC(Cl)=C(F)C(Cl)=C2)C(F)(F)F)C1=O '
               'NAME: (4R,5R)-isocycloseram REASON: CORRECT Contains 4 '
               'fluorine atoms (3 C-F bonds, 3 CF2 groups, 1 CF3 group)\n'
               ' * SMILES: '
               'F[C@@]12C3(C(C(O)CC1C=4CC[C@@](C4CC2O)(O)C)CC(=O)CC3)C NAME: '
               '(1S,9bR)-9b-fluoro-1,5,10-trihydroxy-1,9a-dimethyl-2H,3H,3bH,4H,5H,5aH,6H,8H,9H,10H,11H-cyclopenta[a]phenanthren-7-one '
               'REASON: CORRECT Contains 1 fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               'CC(=O)OCC(=O)[C@@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]12C '
               'NAME: fludrocortisone acetate REASON: CORRECT Contains 1 '
               'fluorine atom (1 C-F bond)\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(N=CC(=C2)C#CC3=CC=C(C=C3)C(F)(F)F)O[C@H]1CN(C)S(=O)(=O)C4=CC=CS4)[C@H](C)CO '
               'NAME: '
               'N-[[(2R,3S)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-8-[2-[4-(trifluoromethyl)phenyl]ethynyl]-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-2-yl]methyl]-N-methyl-2-thiophenesulfonamide '
               'REASON: CORRECT Contains 3 fluorine atoms (3 C-F bonds, 3 CF2 '
               'groups, 1 CF3 group)\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CNC1=NC2=NC(C)=C(C=C2C=N1)C1=CC(NC(=O)NCCC(C)(C)C)=C(F)C=C1C '
               'NAME: LY3009120 REASON: MISSED Contains C and F atoms but no '
               'C-F bonds\n'
               ' * SMILES: '
               'CC(C)n1c(\\C=C\\[C@@H](O)C[C@H](O)CC(O)=O)c(-c2ccc(F)cc2)c2ccccc12 '
               'NAME: '
               '(3S,5S,6E)-7-[3-(4-fluorophenyl)-1-(propan-2-yl)-1H-indol-2-yl]-3,5-dihydroxyhept-6-enoic '
               'acid REASON: MISSED Contains C and F atoms but no C-F bonds\n'
               ' * SMILES: CC(C)N(C(=O)CS(O)(=O)=O)c1ccc(F)cc1 NAME: '
               'flufenacet ESA REASON: MISSED Contains C and F atoms but no '
               'C-F bonds\n'
               ' * SMILES: ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1 '
               'NAME: '
               '3-(4-chlorophenyl)-5-(3-phenyl-1H-pyrazol-4-yl)-1,2,4-oxadiazole '
               'REASON: MISSED Contains no fluorine atoms\n'
               ' * SMILES: '
               'CCC(=O)N1C[C@H](COC[C@@H]2[C@@H]1CC[C@H](O2)CC(=O)NCC3=CC(=CC(=C3)F)F)O '
               'NAME: '
               '2-[(3R,6aS,8S,10aS)-3-hydroxy-1-(1-oxopropyl)-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocin-8-yl]-N-[(3,5-difluorophenyl)methyl]acetamide '
               'REASON: MISSED Contains C and F atoms but no C-F bonds\n'
               ' * SMILES: C=1C(=CC=C(C1)F)C2N(C(CS2)=O)C3=NC=CC=C3 NAME: '
               '2-(4-Fluorophenyl)-3-(2-pyridyl)thiazolidin-4-one REASON: '
               'MISSED Contains C and F atoms but no C-F bonds\n'
               ' * SMILES: '
               'OCC1OC(C(O)C(O)C1O)N1C=C(CC2=CC=C(C=C2)C2CC2)C2=C(F)C=CC=C12 '
               'NAME: TA-1887 REASON: MISSED Contains C and F atoms but no C-F '
               'bonds\n'
               ' * SMILES: COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCCN NAME: '
               'linkable vandetanib analogue REASON: MISSED Contains C and F '
               'atoms but no C-F bonds\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(N=CC(=C2)C#CC3=CC=CC=C3F)O[C@H]1CN(C)C(=O)C4CCCCC4)[C@@H](C)CO '
               'NAME: '
               'N-[[(2R,3R)-8-[2-(2-fluorophenyl)ethynyl]-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-2-yl]methyl]-N-methylcyclohexanecarboxamide '
               'REASON: MISSED Contains C and F atoms but no C-F bonds\n'
               ' * SMILES: '
               'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=CC=C3F)O[C@@H]1CNC)[C@H](C)CO '
               'NAME: '
               '(2R)-2-[(4S,5S)-8-[2-(2-fluorophenyl)ethynyl]-4-methyl-5-(methylaminomethyl)-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-2-yl]-1-propanol '
               'REASON: MISSED Contains C and F atoms but no C-F bonds\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'Contains no fluorine atoms'},
                                 {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(CCC=2C(CC(O)CC2C)(C)C)C)COC3OCC(O)(C3O)CO',
                                     'name': '(3b,9R)-5-Megastigmene-3,9-diol '
                                             '9-[apiosyl-(1->6)-glucoside]',
                                     'reason': 'Contains no fluorine atoms'}],
    'sample_false_negatives': [   {   'smiles': 'ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1',
                                      'name': '3-(4-chlorophenyl)-5-(3-phenyl-1H-pyrazol-4-yl)-1,2,4-oxadiazole',
                                      'reason': 'Contains no fluorine atoms'},
                                  {   'smiles': 'CC1=CC(C)=C(NC(=O)C2=CC(=CC(=C2)C#N)S(F)(F)(F)(F)F)C=C1N1C=CN2N=C(C=C12)C1=CC=CN=C1',
                                      'name': 'BAY-826',
                                      'reason': 'No carbon-fluorine bonds '
                                                'found'},
                                  {   'smiles': 'Nc1cccc2[C@H]3CC[C@H](C3=C(Cl)Cl)c12',
                                      'name': '(1R,4S)-benzovindiflupyr',
                                      'reason': 'Contains no fluorine atoms'}],
    'attempt': 1,
    'reasoning': 'Looking at the outcomes, I notice that the main issue is '
                 'with false negatives - the program is failing to detect some '
                 'valid organofluorine compounds. The key patterns I observe '
                 'in the false negatives are:\n'
                 '\n'
                 '1. The program seems to be missing aromatic C-F bonds, '
                 'particularly in cases like "c1ccc(F)cc1" (fluorobenzene '
                 'patterns)\n'
                 '2. The SMARTS pattern "[C]-[F]" is too specific - it only '
                 'matches aliphatic C-F single bonds\n'
                 '3. Some compounds have fluorine atoms that are not being '
                 'properly detected by the current SMARTS patterns\n'
                 '\n'
                 "Here's how I'll improve the program:\n"
                 '\n'
                 '1. Add aromatic carbon-fluorine bond pattern with "c-F" '
                 'SMARTS\n'
                 '2. Make the C-F bond pattern more general by using "C,c" to '
                 'match both aliphatic and aromatic carbons\n'
                 '3. Simplify the approach to just look for any '
                 'carbon-fluorine bonds rather than trying to distinguish '
                 'between CF, CF2 and CF3 groups\n'
                 '4. Add additional validation to ensure the detected bonds '
                 'are actually C-F bonds\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1229,
    'num_false_positives': 100,
    'num_true_negatives': 1900,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.9247554552294959,
    'recall': 0.997564935064935,
    'f1': 0.9597813354158531,
    'accuracy': 0.9681311881188119,
    'negative_predictive_value': 0.9984235417761429}