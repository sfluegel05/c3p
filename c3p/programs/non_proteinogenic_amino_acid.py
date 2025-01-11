"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determine if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    A non-proteinogenic amino acid is defined as one with both amino and carboxyl groups,
    plus non-standard modifications in side chains or additional groups not seen in
    the 20 standard amino acids.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for primary, secondary, and tertiary amines
    amino_group_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0]')
    # SMARTS pattern for carboxyl group
    carboxyl_group_pattern = Chem.MolFromSmarts('C(=O)[O-,O]')
    
    # Check for amino and carboxyl groups presence
    has_amino_group = mol.HasSubstructMatch(amino_group_pattern)
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_group_pattern)
    
    if not (has_amino_group and has_carboxyl_group):
        return False, "Must contain both amino and carboxyl groups"

    # SMARTS for oxidative modifications and added complexity in side chains
    unique_modifications = [
        Chem.MolFromSmarts('C=C'),    # Double bonds in unusual places
        Chem.MolFromSmarts('N=C'),    # Amide or imine modifications
        Chem.MolFromSmarts('N-[NX3]'), # Diverse nitrogen reduction functionality
        Chem.MolFromSmarts('S'),      # Sulfur atoms outside standard placement
        Chem.MolFromSmarts('[F,Cl,Br,I]'), # Halogenation patterns
        Chem.MolFromSmarts('O=C-N'),  # Acetyl groups or modified amidation
        Chem.MolFromSmarts('[n,N]'),  # Aromatic nitrogens
        Chem.MolFromSmarts('[CH2]O[CH2]'), # Ether linkages
        Chem.MolFromSmarts('CCC(=O)'), # Keto groups in side chains
        Chem.MolFromSmarts('c1nccc1'), # Uncommon aromatic rings
    ]

    # Check for these unique side chain patterns
    for pattern in unique_modifications:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unique/modified side chain structures indicating non-proteinogenic nature"

    return False, "Does not have unique/modified side chain features expected in non-proteinogenic amino acids"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83820',
                          'name': 'non-proteinogenic amino acid',
                          'definition': 'Any amino-acid that is not naturally '
                                        'encoded in the genetic code of any '
                                        'organism.',
                          'parents': ['CHEBI:33709'],
                          'xrefs': ['Wikipedia:Non-proteinogenic_amino_acids'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7222222222222222 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C(/C=C\\CCCCC)[C@H](C/C=C/C=C/[C@H]([C@@H](O)CCCC(=O)O)SC[C@H](N)C(=O)O)O '
               'NAME: (12R)-hydroxy-10,11-dihydroleukotriene E4 REASON: '
               'CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               ' * SMILES: OC([C@H](CCS(O)(=O)=O)N)=O NAME: L-homocysteic acid '
               'REASON: CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               ' * SMILES: N[14C@]1(C[C@H](F)C1)C(O)=O NAME: fluciclovine '
               '((14)C) REASON: CORRECT Contains unique/modified side chain '
               'structures indicating non-proteinogenic nature\n'
               ' * SMILES: Cc1nc(NCCC[C@H](N)C(O)=O)nc(C)c1O NAME: '
               'argpyrimidine REASON: CORRECT Contains unique/modified side '
               'chain structures indicating non-proteinogenic nature\n'
               ' * SMILES: O=C(O)[C@@H](NO)CCCCSC NAME: '
               'N-hydroxy-L-dihomomethionine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: NC(CCCNC(N)=O)C(O)=O NAME: citrulline REASON: '
               'CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               ' * SMILES: N\\C(=C/C=C/C(O)=O)C(O)=O NAME: '
               '(2Z,4E)-2-aminomuconic acid REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: C(N(C(=NO)NC)O)CC[C@@H](C(=O)[O-])[NH3+] NAME: '
               "N(delta),N(omega')-dihydroxy-N(omega)-methyl-L-arginine "
               'zwitterion REASON: CORRECT Contains unique/modified side chain '
               'structures indicating non-proteinogenic nature\n'
               ' * SMILES: '
               'C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)[C@H](N)C(O)=O '
               'NAME: FMN-L-threonine REASON: CORRECT Contains unique/modified '
               'side chain structures indicating non-proteinogenic nature\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](NC(=O)C[C@H](N)C(O)=O)[C@H](NC(C)=O)[C@@H](O)[C@@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O '
               'NAME: '
               'N(4)-{N-acetyl-beta-D-glucosaminyl-(1->2)-alpha-D-mannosyl-(1->3)-[N-acetyl-beta-D-glucosaminyl-(1->2)-alpha-D-mannosyl-(1->6)]-beta-D-mannosyl-(1->4)-N-acetyl-beta-D-glucosaminyl-(1->4)-[alpha-L-fucosyl-(1->6)]-N-acetyl-beta-D-glucosaminyl}-L-asparagine '
               'REASON: CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               ' * SMILES: N[C@@H](Cn1ccc(=O)o1)C(O)=O NAME: '
               '3-(5-oxoisoxazolin-2-yl)-L-alanine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: CSCCCCCC(N(O)O)C(O)=O NAME: '
               'N,N-dihydroxytrihomomethionine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: '
               'N[C@@H](CCSC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(O)ncnc12)C(O)=O '
               'NAME: S-inosyl-L-homocysteine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: O=C([O-])[C@@H]([NH3+])[C@@H](CCNC(=[NH2+])N)C '
               'NAME: (3R)-3-methyl-L-arginine zwitterion REASON: CORRECT '
               'Contains unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: CSC(CC(N)C(O)=O)C(O)=O NAME: gamma-thiomethyl '
               'glutamate REASON: CORRECT Contains unique/modified side chain '
               'structures indicating non-proteinogenic nature\n'
               ' * SMILES: CC(C)=CCSC[C@H](N)C(O)=O NAME: S-prenyl-L-cysteine '
               'REASON: CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               ' * SMILES: C(C(C(O)=O)N)C1=CC=CC=C1F NAME: '
               '2-fluorophenylalanine REASON: CORRECT Contains unique/modified '
               'side chain structures indicating non-proteinogenic nature\n'
               ' * SMILES: '
               'CC1(C)C[C@@H]1C(=O)N\\C(C(O)=O)=C/CCCCSC[C@H](N)C(O)=O NAME: '
               'cilastatin REASON: CORRECT Contains unique/modified side chain '
               'structures indicating non-proteinogenic nature\n'
               ' * SMILES: N[C@@H](CCCNC(=N)NO)C(O)=O NAME: '
               'N(5)-[(hydroxyamino)(imino)methyl]-L-ornithine REASON: CORRECT '
               'Contains unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: S(CC[C@@](N)(C)C(O)=O)C NAME: '
               'alpha-methylmethionine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: C1=C(CC(C=C1)N)C(O)=O NAME: '
               '5-aminocyclohexa-1,3-diene-1-carboxylic acid REASON: CORRECT '
               'Contains unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: N[C@@H](CCSC(F)(F)F)C(O)=O NAME: '
               'trifluoro-L-methionine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: CNC(=O)C[C@H](N)C(O)=O NAME: '
               'N(4)-methyl-L-asparagine REASON: CORRECT Contains '
               'unique/modified side chain structures indicating '
               'non-proteinogenic nature\n'
               ' * SMILES: O=C(O)[C@@H](N)CCCCCCCCSC NAME: '
               'L-hexahomomethionine REASON: CORRECT Contains unique/modified '
               'side chain structures indicating non-proteinogenic nature\n'
               ' * SMILES: '
               'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1NC(=O)C[C@H](N)C(O)=O '
               'NAME: N(4)-(beta-N-acetyl-D-glucosaminyl)-L-asparagine REASON: '
               'CORRECT Contains unique/modified side chain structures '
               'indicating non-proteinogenic nature\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O(C1=C(C(=C(C[C@@H](N)C(O)=O)C=C1)C)C)C NAME: '
               'D-O,2,3-trimethyltyrosine REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               ' * SMILES: N[C@H](Cc1ccc(O)c(O)c1)C(O)=O NAME: D-dopa REASON: '
               'MISSED Does not have unique/modified side chain features '
               'expected in non-proteinogenic amino acids\n'
               ' * SMILES: C[C@H](N(C)C)C(O)=O NAME: N,N-dimethyl-L-alanine '
               'REASON: MISSED Does not have unique/modified side chain '
               'features expected in non-proteinogenic amino acids\n'
               ' * SMILES: N[C@@H](C(O)CC(O)=O)C(O)=O NAME: '
               '3-hydroxy-L-glutamic acid REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               ' * SMILES: OCC(Nc1ccc(cc1[N+]([O-])=O)[N+]([O-])=O)C(O)=O '
               'NAME: N-(2,4-dinitrophenyl)serine REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               ' * SMILES: CCCC[C@H](N)C(O)=O NAME: L-norleucine REASON: '
               'MISSED Does not have unique/modified side chain features '
               'expected in non-proteinogenic amino acids\n'
               ' * SMILES: C1=C(C(=C(C=C1C[C@@H](C(O)=O)N)O)OC)C NAME: '
               '3-hydroxy-O,5-dimethyl-L-tyrosine REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               ' * SMILES: CCCCCC(N)C(O)=O NAME: 2-aminoheptanoic acid REASON: '
               'MISSED Does not have unique/modified side chain features '
               'expected in non-proteinogenic amino acids\n'
               ' * SMILES: N[C@@H](CCOP(O)(O)=O)C(O)=O NAME: '
               'O-phospho-L-homoserine REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               ' * SMILES: C1(C(C[C@@H](C(=O)O)N)=O)=CC=CC(=C1N)O NAME: '
               '3-hydroxy-L-kynurenine REASON: MISSED Does not have '
               'unique/modified side chain features expected in '
               'non-proteinogenic amino acids\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'Must contain both amino and '
                                               'carboxyl groups'}],
    'sample_false_negatives': [   {   'smiles': 'C[N+]1(C)CCC[C@H]1C([O-])=O',
                                      'name': 'L-proline betaine',
                                      'reason': 'Must contain both amino and '
                                                'carboxyl groups'},
                                  {   'smiles': 'C[n+]1ccccc1C([O-])=O',
                                      'name': 'Betaine homarine',
                                      'reason': 'Must contain both amino and '
                                                'carboxyl groups'},
                                  {   'smiles': 'C(C(=O)CSSC[C@@H](C(=O)[O-])[NH3+])(=O)[O-]',
                                      'name': '(2R)-2-ammonio-2-carboxylatoethyl-disulfanyl-oxopropanoate',
                                      'reason': 'Must contain both amino and '
                                                'carboxyl groups'},
                                  {   'smiles': 'C([C@@H](C([O-])=O)[NH3+])C=1N(C=NC1)C',
                                      'name': 'N(pros)-methyl-L-histidine '
                                              'zwitterion',
                                      'reason': 'Must contain both amino and '
                                                'carboxyl groups'},
                                  {   'smiles': 'CC(C)(S)[C@H](C(O)=O)n1ccnc1Cc1ccccc1',
                                      'name': 'benzylpenillamine',
                                      'reason': 'Must contain both amino and '
                                                'carboxyl groups'}],
    'attempt': 2,
    'reasoning': 'To improve the function for classifying non-proteinogenic '
                 'amino acids, we need to address the shortcomings of the '
                 'previous implementation, particularly the missed false '
                 "negatives. Here's a step-by-step approach to refining the "
                 'logic and criteria used for classification:\n'
                 '\n'
                 '1. **Revisiting Basic Functional Groups**: Ensure that the '
                 'detection of amino and carboxyl groups is robust. The amino '
                 'group detection should allow more flexibility to capture '
                 'secondary and tertiary amines frequently found in '
                 'non-standard amino acids.\n'
                 '\n'
                 '2. **Diverse Modifications Inclusion**: Increase the variety '
                 'of modifications in side chains considered unique or '
                 'non-standard. This can include broadening the scope for '
                 'halogens, sulfur atoms (beyond known cysteine and methionine '
                 'patterns), ketones, aldehydes, unusual unsaturation '
                 'patterns, large aliphatic chains, or multiple carboxyl or '
                 'amino groups.\n'
                 '\n'
                 '3. **Mismatch Correction**: Ensure the molecule is '
                 'categorized correctly if it matches other criteria but lacks '
                 'clear unconventional fingerings. There should be provisions '
                 'to better capture subtle modifications, such as ring '
                 'oxidation or extended conjugation chains.\n'
                 '\n'
                 '4. **Reactivity and Environment Recognition**: Reflect '
                 'changes and reactivity that differentiate non-proteinogenic '
                 'amino acids from the standard set, like ringed structures '
                 'that are not in standard amino acids, or additional nitrogen '
                 'functionality suggesting cyclization or novel biochemistry '
                 'roles.\n'
                 '\n'
                 "Here's the updated code to achieve better classification:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 678,
    'num_false_positives': 100,
    'num_true_negatives': 441,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.87146529562982,
    'recall': 0.9926793557833089,
    'f1': 0.9281314168377824,
    'accuracy': 0.9142156862745098,
    'negative_predictive_value': 0.9887892376681614}