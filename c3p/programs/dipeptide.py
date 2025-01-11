"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is likely to be a dipeptide based on its SMILES string.
    A dipeptide is generally composed of two amino acid units, possibly within more complex connections.
    
    It is vital to recognize varying structures, allowing some flexibility in peptide bond counts,
    provided there are recognizably two amino acid-like structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide variant, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for identifying peptide bonds (-CO-NH- linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Allow for one or two peptide bonds to account for variability in structure
    if len(peptide_bond_matches) < 1:
        return False, f"Expected at least 1 peptide bond, found {len(peptide_bond_matches)}"

    # SMARTS pattern for a more generalized amino acid core
    amino_acid_core_pattern = Chem.MolFromSmarts("[NX3][CX4]C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_core_pattern)
    
    # Check for at least two amino acid-like structures
    if len(amino_acid_matches) < 2:
        return False, f"Expected at least 2 amino acid residues, found {len(amino_acid_matches)}"

    return True, "Molecule likely contains two amino acid units that form a dipeptide structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46761',
                          'name': 'dipeptide',
                          'definition': 'Any molecule that contains two '
                                        'amino-acid residues connected by '
                                        'peptide linkages.',
                          'parents': ['CHEBI:25676'],
                          'xrefs': ['KEGG:C00107'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               '[H][C@@]12CCCN1C(=O)[C@H](Cc1c([nH]c3ccccc13)C(C)(C)C=C)NC2=O '
               'NAME: deoxybrevianamide E REASON: MISSED Expected 1 peptide '
               'bond, found 2\n'
               ' * SMILES: C[C@@H](N)C(=O)N[C@H](C)C(O)=O NAME: '
               'D-alanyl-D-alanine REASON: MISSED Expected 2 amino acid '
               'residues, found 1\n'
               ' * SMILES: NCC(=O)N[C@@H](CCCNC(N)=N)C(O)=O NAME: Gly-Arg '
               'REASON: MISSED Expected 2 amino acid residues, found 1\n'
               ' * SMILES: C([C@H]1NC(CC1)=O)(=O)NCC(O)=O NAME: '
               'pyroglutamylglycine REASON: MISSED Expected 1 peptide bond, '
               'found 2\n'
               ' * SMILES: '
               'C(N[C@@H](C(C)C)C(=O)N[C@@H](CC(C)C)C(N[C@H](/C=C/C(OCC)=O)C[C@H]1C(NCC1)=O)=O)(OCC2=CC=CC=C2)=O '
               'NAME: Ethyl '
               '(E,4S)-4-[[(2S)-4-methyl-2-[[(2S)-3-methyl-2-(phenylmethoxycarbonylamino)butanoyl]amino]pentanoyl]amino]-5-[(3S)-2-oxopyrrolidin-3-yl]pent-2-enoate '
               'REASON: MISSED Expected 1 peptide bond, found 4\n'
               ' * SMILES: '
               'O=C(N1C(CCC1)C(=O)NC2=C(C=CC=C2O)C(O)=O)C(NC(=O)C3=CC=C(O)C=C3)C(C)C '
               'NAME: '
               '3-hydroxy-2-[[1-[2-[(4-hydroxybenzoyl)amino]-3-methylbutanoyl]pyrrolidine-2-carbonyl]amino]benzoic '
               'acid REASON: MISSED Expected 1 peptide bond, found 3\n'
               ' * SMILES: '
               'O=C1NCC[C@H](O)C[C@H](NC(=O)[C@@H](NC(=O)/C=C/C=C/CCCCCCCCC)[C@H](O)C)C(N[C@H](C=C1)C)=O '
               'NAME: Glidobactin C REASON: MISSED Expected 1 peptide bond, '
               'found 4\n'
               ' * SMILES: N[C@@H](Cc1c[nH]cn1)C(=O)NCC(O)=O NAME: His-Gly '
               'REASON: MISSED Expected 2 amino acid residues, found 1\n'
               ' * SMILES: N[C@H](CCC(=O)N[C@H](CCC(O)=O)C(O)=O)C(O)=O NAME: '
               'D-gamma-glutamyl-D-glutamic acid REASON: MISSED Expected 2 '
               'amino acid residues, found 1\n'
               ' * SMILES: '
               'S1SC23N(C4C(OC(=O)C)C=COC=C4C2)C(C15N(C6C(OC(=O)C)C=COC=C6C5)C3=O)=O '
               'NAME: Acetylaranotin REASON: MISSED Expected 1 peptide bond, '
               'found 2\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Expected at least 1 peptide '
                                               'bond, found 0'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'Expected at least 2 amino acid '
                                               'residues, found 0'},
                                 {   'smiles': 'C1C[C@@H]([C@@H](O[C@@H]1CCNC(=O)C2CCOCC2)CO)NC(=O)NC3=CC(=CC=C3)Cl',
                                     'name': 'N-[2-[(2S,5S,6R)-5-[[(3-chloroanilino)-oxomethyl]amino]-6-(hydroxymethyl)-2-oxanyl]ethyl]-4-oxanecarboxamide',
                                     'reason': 'Expected at least 2 amino acid '
                                               'residues, found 0'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)[C@H]([NH2+]C)C(=O)N[C@@H](CC=1C=2C=CC=CC2NC1)CO',
                                      'name': 'N-methyl-L-valyl-L-tryptophanol(1+)',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 0'},
                                  {   'smiles': 'O=C\\1N[C@H](C(=O)N/C1=C/C=2C=3C(NC2C(C=C)(C)C)=CC=C(C3)CC=C(C)C)C',
                                      'name': 'isoechinulin A',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'C[C@H](NC(=O)CCN)C(O)=O',
                                      'name': 'beta-Ala-Ala',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'O=C\\1NC(C(=O)N/C1=C/C=2C=3C(NC2C(C=C)(C)C)=CC=C(C3)CC=C(C)C)=C',
                                      'name': 'isoechinulin B',
                                      'reason': 'Expected at least 1 peptide '
                                                'bond, found 0'},
                                  {   'smiles': 'O=C(NCCC(=O)O)[C@@H](N(C1=CC=C(C=2C1=NON2)[N+]([O-])=O)C)CC(=O)O',
                                      'name': 'NBD-N(Me)-DbetaA-OH',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'O=C\\1NC(C(=O)[N-]/C1=C/C=2C=3C(NC2C(C=C)(C)C)=CC=C(C3)CC=C(C)C)=C',
                                      'name': 'isoechinulin B anion',
                                      'reason': 'Expected at least 1 peptide '
                                                'bond, found 0'},
                                  {   'smiles': 'N[C@@H](Cc1ccccc1)C(=O)NCCC(O)=O',
                                      'name': 'Phe-beta-Ala',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'CC(=O)NCCC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O',
                                      'name': 'N-acetylcarnosine',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'OC(=O)c1ccc(NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)c2ccccc2)cc1',
                                      'name': 'bentiromide',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'},
                                  {   'smiles': 'O=C(NC(C(CC)C)C([O-])=O)C([NH3+])C(C)C',
                                      'name': '2-[(2-Azaniumyl-3-methylbutanoyl)amino]-3-methylpentanoate',
                                      'reason': 'Expected at least 2 amino '
                                                'acid residues, found 1'}],
    'attempt': 2,
    'reasoning': 'To improve the classification of dipeptides based on the '
                 "SMILES string, let's first analyze the previous program and "
                 'outcomes:\n'
                 '\n'
                 '1. **Issue**: The expectation was a single peptide bond '
                 '(`-C(=O)N-` linkage pattern), but some molecules have '
                 'multiple peptide bonds.\n'
                 '   - **Solution**: Adjust the checker to allow for two '
                 'peptide bonds, which can occur in cyclic peptides or when '
                 'two residues are connected via additional peptide bonds.\n'
                 '   \n'
                 '2. **Issue**: The program expects exactly two amino acid '
                 'residues. However, the peptide may sometimes have more '
                 'complex or cyclic arrangements, needing a more flexible '
                 'approach.\n'
                 '   - **Solution**: Check for at least two amino acid '
                 'patterns and ensure the overall molecular structure aligns '
                 'with a dipeptide concept.\n'
                 '\n'
                 '3. **Issue**: The SMARTS patterns for amino acids were too '
                 'specific, missing correctly arranged peptides.\n'
                 '   - **Solution**: Use a more general pattern to capture '
                 'amino acid-like structures while confirming at least two '
                 'such structures exist.\n'
                 '\n'
                 "Let's adjust the code considering these insights:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 924,
    'num_false_positives': 100,
    'num_true_negatives': 1052,
    'num_false_negatives': 45,
    'num_negatives': None,
    'precision': 0.90234375,
    'recall': 0.9535603715170279,
    'f1': 0.9272453587556447,
    'accuracy': 0.9316360207449317,
    'negative_predictive_value': 0.95897903372835}