"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol typically contains the -SH group connected to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a thiol group (covering open chains, rings, etc.)
    thiol_patterns = [
        Chem.MolFromSmarts("[SX2H]"),          # standard thiol
        Chem.MolFromSmarts("[cH][SX2H]"),      # thiol adjacent to aromatic carbon
        Chem.MolFromSmarts("[C;!R]S"),         # aliphatic thiol
    ]

    # Look for thiol group in the molecule
    for pattern in thiol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a thiol group (-SH)"

    return False, "No thiol group found"


# Example usage
smiles_examples = [
    "CCSCCS", "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCS", "SCC(CC)C",
    "CCO[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CCS(C)(=O)=O)C(=O)OC(C)C",
    "NC(=O)CCCCC(S)CCS", "SC(CS)(C)C", "Mary had a little lamb",
    "SCCCCCCCS", "NCCS"
]

# Test the function
for sml in smiles_examples:
    result, reason = is_thiol(sml)
    print(f"SMILES: {sml} => is_thiol: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29256',
                          'name': 'thiol',
                          'definition': 'An organosulfur compound in which a '
                                        'thiol group, -SH, is attached to a '
                                        'carbon atom of any aliphatic or '
                                        'aromatic moiety.',
                          'parents': ['CHEBI:33261'],
                          'xrefs': ['KEGG:C00145', 'Wikipedia:Thiol'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5492957746478873 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CCSCCS NAME: '
               '2-(ethylsulfanyl)ethanethiol REASON: CORRECT Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCS NAME: '
               'pantetheine REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SCC(CC)C NAME: 2-Methyl-1-butanethiol REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: '
               'CC[C@@H](C)[C@@H](CO[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CCS(C)(=O)=O)C(=O)OC(C)C)NCC(CS)N '
               'NAME: '
               '(2S)-2-[[(2S)-2-[(2S,3R)-2-[(2-amino-3-mercaptopropyl)amino]-3-methylpentoxy]-1-oxo-3-phenylpropyl]amino]-4-methylsulfonylbutanoic '
               'acid propan-2-yl ester REASON: CORRECT Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: NC(=O)CCCCC(S)CCS NAME: dihydrolipoamide REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SC(CS)(C)C NAME: 2-methylpropane-1,2-dithiol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SCCCC(O)=O NAME: 4-sulfanylbutanoic acid REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SC(C(CO)C)CC NAME: 3-Mercapto-2-methylpentanol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: NCCS NAME: cysteamine REASON: CORRECT Contains a '
               'thiol group (-SH)\n'
               ' * SMILES: SCCCCCCCCCS NAME: 1,9-Nonanedithiol REASON: CORRECT '
               'Contains a thiol group (-SH)\n'
               ' * SMILES: SC(C(C)C)CO NAME: xi-2-Mercapto-3-methyl-1-butanol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: NCCCNCCS NAME: WR-1065 REASON: CORRECT Contains a '
               'thiol group (-SH)\n'
               ' * SMILES: SC/C=C(/CCC=C(C)C)\\C NAME: thiogeraniol REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SCCCS NAME: propane-1,3-dithiol REASON: CORRECT '
               'Contains a thiol group (-SH)\n'
               ' * SMILES: OC(=O)C(=O)CS NAME: 3-mercaptopyruvic acid REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: C=1C=CC=C2NC=C(C12)CCCC(N[C@H](C(=O)O)CS)=O NAME: '
               'N-[4-(indol-3-yl)butanoyl]-L-cysteine REASON: CORRECT Contains '
               'a thiol group (-SH)\n'
               ' * SMILES: CC(C)(S)[C@H](N)C(O)=O NAME: L-penicillamine '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: CC(C)(S)C(N)C(O)=O NAME: penicillamine REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: '
               'CC[C@H](C)[C@@H](COC(CC1=CC=CC=C1)C(=O)N[C@@H](CCS(C)(=O)=O)C(=O)OC(C)C)NC[C@@H](N)CS '
               'NAME: L-744,832 REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: OC(CS)C(O)CS NAME: 1,4-dimercaptobutane-2,3-diol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SC(CCC)(CO)C NAME: 2-Mercapto-2-methyl-1-pentanol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SCC=C(C)C NAME: 3-Methyl-2-butene-1-thiol REASON: '
               'CORRECT Contains a thiol group (-SH)\n'
               ' * SMILES: SC(CCCCC)C NAME: 2-Heptanethiol REASON: CORRECT '
               'Contains a thiol group (-SH)\n'
               ' * SMILES: OC(=O)C(CS)CS NAME: '
               '3-mercapto-2-mercaptomethylpropanoic acid REASON: CORRECT '
               'Contains a thiol group (-SH)\n'
               ' * SMILES: O[C@H](CS)[C@H](O)CS NAME: D-1,4-dithiothreitol '
               'REASON: CORRECT Contains a thiol group (-SH)\n'
               'False positives: SMILES: '
               'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO NAME: '
               'Cys-Ser-Gln REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CCC(=O)N)C(=O)N[C@@H](CCCCN)C(O)=O '
               'NAME: Gln-Cys-Lys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: SC[C@H](NC(=O)CNC(=O)[C@@H](N)[C@H](O)C)C(O)=O '
               'NAME: Thr-Gly-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: SCC(NC(=O)C(N)CC(O)=O)C(O)=O NAME: '
               'Aspartyl-Cysteine REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CO)C(=O)N[C@@H](CC(=O)N)C(O)=O NAME: '
               'Ser-Cys-Asn REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: SC[C@H](NC(=O)[C@@H](NC(=O)[C@H]1NCCC1)CCCCN)C(O)=O '
               'NAME: Pro-Lys-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CCC(=O)N)CCC(=O)N)C(O)=O '
               'NAME: Gln-Gln-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'CNC(=O)C(NC(=O)[C@H](CC(C)C)NC(=O)C(S)CCN1C(=O)CCC1=O)C(C)(C)C '
               'NAME: D1927 REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(=O)N[C@@H](CCC(=O)N)C(O)=O '
               'NAME: Cys-Trp-Gln REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'N[C@@H](CCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(O)=O)C(O)=O)C(O)=O)C(O)=O '
               'NAME: phytochelatin 3 REASON: WRONGLY CLASSIFIED Contains a '
               'thiol group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCCN=C(N)N)C(O)=O)CC(O)=O '
               'NAME: Cys-Asp-Arg REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@@H]([C@H](O)C)C(=O)N[C@@H](CCC(=O)N)C(O)=O '
               'NAME: Cys-Thr-Gln REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(O)=O '
               'NAME: Cys-Arg-Tyr REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N1[C@@H](CCC1)C(=O)N[C@@H](CC2=CC=CC=C2)C(O)=O '
               'NAME: Cys-Pro-Phe REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CCC(O)=O)CCC(=O)N)C(O)=O '
               'NAME: Glu-Gln-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CS)C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(O)=O '
               'NAME: Cys-Cys-Tyr REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CCSC)C(=O)N[C@@H](CCC(O)=O)C(O)=O NAME: '
               'Met-Cys-Glu REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C)CCC(O)=O)C(O)=O NAME: '
               'Ala-Glu-Cys REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: '
               'ClC1=C2N(C(=O)C[C@H](OC(=O)[C@@H](N(C)C(=O)CCC(S)(C)C)C)[C@@]3(O[C@H]3[C@@H]([C@H]4OC(=O)N[C@](O)(C4)[C@H](OC)C=CC=C(CC(=C2)C=C1OC)C)C)C)C '
               'NAME: Ravtansine REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CCCCN)C(O)=O '
               'NAME: Arg-Lys-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CCC(O)=O)C(O)=O '
               'NAME: Tyr-Glu-Cys REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CS)N[C@H](C(N1[C@H](C(=O)N[C@H](C(O)=O)CC=2C=CC(=CC2)O)CCC1)=O)CCC(=O)N '
               'NAME: Cys-Gln-Pro-Tyr REASON: WRONGLY CLASSIFIED Contains a '
               'thiol group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](N)C(=O)N[C@@H](CC1=CC=C(O)C=C1)C(=O)N[C@@H](C(C)C)C(O)=O '
               'NAME: Cys-Tyr-Val REASON: WRONGLY CLASSIFIED Contains a thiol '
               'group (-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC(=O)N)CCCCN)C(O)=O NAME: '
               'Asn-Lys-Cys REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               ' * SMILES: '
               'SC[C@H](NC(=O)[C@@H](N)CS)C(=O)N[C@@H](CC(C)C)C(O)=O NAME: '
               'Cys-Cys-Leu REASON: WRONGLY CLASSIFIED Contains a thiol group '
               '(-SH)\n'
               'False negatives: SMILES: '
               'CC(C)=CCNc1nc(S)nc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O '
               "NAME: 2-thio-N(6)-dimethylallyladenine 5'-monophosphate "
               'REASON: MISSED No thiol group found\n'
               ' * SMILES: SC=1C(P)=CC=CC1 NAME: 2-Phosphinobenzenethiol '
               'REASON: MISSED No thiol group found\n'
               ' * SMILES: '
               'O[C@](CN1N=CN=C1S)(C1=CC=C(F)C=C1F)C(F)(F)C1=CC=C(OC2=CC=C(C=C2)C#N)C=N1 '
               'NAME: (S)-fluoxytioconazole (thiol) REASON: MISSED No thiol '
               'group found\n'
               ' * SMILES: '
               'C[C@@H](O)[C@H](NC(=O)Nc1nc(S)nc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O)C(O)=O '
               "NAME: 2-thio-N(6)-L-threonylcarbamoyladenine 5'-monophosphate "
               'REASON: MISSED No thiol group found\n'
               ' * SMILES: Sc1ncnc2nc[nH]c12 NAME: purine-6-thiol REASON: '
               'MISSED No thiol group found\n'
               ' * SMILES: SC=1C(CC)=CC=CC1 NAME: 2-Ethylbenzenethiol REASON: '
               'MISSED No thiol group found\n'
               ' * SMILES: OC(=O)c1ccccc1S NAME: thiosalicylic acid REASON: '
               'MISSED No thiol group found\n'
               ' * SMILES: Cn1cnc(S)c1C[C@H](N)C(O)=O NAME: ovothiol A REASON: '
               'MISSED No thiol group found\n'
               ' * SMILES: COC1=C(C=CC(=C1)SC)C(=O)NC2=C(C=C(C=N2)Cl)Cl NAME: '
               'N-(3,5-dichloro-2-pyridinyl)-2-methoxy-4-(methylthio)benzamide '
               'REASON: MISSED No thiol group found\n'
               ' * SMILES: SC1=CC=C(F)C=C1 NAME: 4-Fluorothiophenol REASON: '
               'MISSED No thiol group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No thiol group found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No thiol group found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]12Nc3nc(N)[nH]c(=O)c3N[C@@]1([H])C1=C(S[Mo](O)(O)(O)S1)[C@@H](COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(N)nc1=O)O2',
                                      'name': 'Mo(V)-molybdopterin cytosine '
                                              'dinucleotide',
                                      'reason': 'No thiol group found'}],
    'attempt': 1,
    'reasoning': "To address the previous failure, let's first analyze the "
                 'outcomes:\n'
                 '\n'
                 '- **True Positives:** The function correctly detected many '
                 'thiol-containing molecules. This indicates the pattern for '
                 'thiol detection worked in these cases.\n'
                 '  \n'
                 '- **False Positives:** The test classified some peptides '
                 'with a cysteine (which contains a thiol group) as having a '
                 'thiol group. While these peptides do contain a thiol group, '
                 'the classification appears overly permissive in including '
                 'them perhaps due to incorrect substructure matches.\n'
                 '\n'
                 '- **False Negatives:** The function missed other '
                 'thiol-containing compounds (e.g., some specific aromatic '
                 'thiols and thiol-containing cyclic structures). This '
                 'suggests that the pattern used to identify thiols was too '
                 'restrictive or did not cover all forms in which thiol groups '
                 'can appear within complex structures.\n'
                 '\n'
                 "Here's a revised plan to improve the classifier:\n"
                 '\n'
                 "1. **Optimize SMARTS Pattern:** We'll look for thiol groups "
                 'not just as peripheral atoms but as part of more complex '
                 'motifs. Expanding the pattern to include more variations '
                 'could help reduce false negatives.\n'
                 '\n'
                 '2. **Redefine the pattern for thiol detection**, potentially '
                 'focusing on different hybridization states and geometries '
                 'where sulfur might exist (e.g., identifying aromatic sulfur '
                 'or different valency sulfur in conjugated systems).\n'
                 '\n'
                 '3. **Restrict misclassifications due to peptides:** We need '
                 'better filters or logic to prevent misidentifying certain '
                 'thiol-containing substructures in larger peptide chains as '
                 'standalone thiol groups, which might help reduce false '
                 'positives.\n'
                 '\n'
                 '4. **Consider sulfur as part of different substructures:** '
                 'We need to capture thiols as part of rings or more complex '
                 'frameworks accurately.\n'
                 '\n'
                 'With these improvements in mind, here is the revised code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 105,
    'num_false_positives': 100,
    'num_true_negatives': 1283,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.5121951219512195,
    'recall': 0.9905660377358491,
    'f1': 0.6752411575562701,
    'accuracy': 0.9321692411014103,
    'negative_predictive_value': 0.9992211838006231}