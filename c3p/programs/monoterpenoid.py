"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:35189 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes (C10 skeleton) but may have some
    carbons removed or rearranged, and often contain oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Core structural patterns for monoterpenoid skeleton
    core_patterns = {
        "isopropyl": "[CH3][CH1]([CH3])[CH1,CH2]",
        "gem_dimethyl": "[CH3][C]([CH3])",
        "isoprene_unit": "[CH3,CH2][C]([CH3,CH2])=[C,c]",
        "menthane_skeleton": "CC1CCC(C(C)C)CC1",
        "pinane_skeleton": "C12CCC(C1(C)C)C2",
        "carane_skeleton": "C12CC1C(CC2)",
        "thujane_skeleton": "CC(C)C12CC1CC2",
        "fenchane_skeleton": "CC1(C)C2CCC1(C)C2",
        "bornane_skeleton": "CC1(C)C2CCC1(C)CC2",
    }
    
    # Additional patterns for common modifications
    mod_patterns = {
        "cyclohexane": "C1CCCCC1",
        "cyclopentane": "C1CCCC1",
        "bridged_bicycle": "C12CCC1CC2",
        "aromatic_ring": "c1ccccc1",
        "conjugated_system": "C=CC=C",
    }
    
    # Functional group patterns
    func_patterns = {
        "alcohol": "[OH1]",
        "ether": "[OR0]",
        "ketone": "[CX3](=[OX1])[#6]",
        "ester": "[#6][CX3](=[OX1])[OX2][#6]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "epoxide": "C1OC1",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "enol": "[OX2H1][#6X3]=[#6X3]",
        "hemiacetal": "[OX2H1][CX4]([OX2])",
        "amino": "[NX3;H2,H1;!$(NC=O)]",
        "amide": "[NX3][CX3](=[OX1])[#6]",
    }
    
    # Check for core patterns
    core_features = []
    for name, smarts in core_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            core_features.append(name.replace("_", " "))

    # Check for modification patterns
    mod_features = []
    for name, smarts in mod_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            mod_features.append(name.replace("_", " "))

    # Check for functional groups
    func_features = []
    for name, smarts in func_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            func_features.append(name)

    # Ring analysis
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(r) for r in ring_info.AtomRings()]
    
    # Score calculation
    score = 0
    
    # Core features are most important
    score += len(core_features) * 2
    
    # Modification patterns add support
    score += len(mod_features)
    
    # Functional groups are common in monoterpenoids
    if func_features:
        score += min(len(func_features), 3)  # Cap at 3 to avoid over-counting
    
    # Ring count consideration
    if 1 <= ring_count <= 4:
        score += 1
    
    # Common ring sizes (5,6) are typical
    if any(rs in [5,6] for rs in ring_sizes):
        score += 1
    
    # Molecular weight consideration (more flexible range)
    if 100 <= mol_wt <= 500:
        score += 1

    # Combine all features for reporting
    all_features = []
    if core_features:
        all_features.extend(core_features)
    if mod_features:
        all_features.extend(mod_features)
    if func_features:
        all_features.append("oxygen/nitrogen-containing groups")
    if ring_count > 0:
        all_features.append(f"{ring_count} ring(s)")

    # Classification
    # Require either strong core features or combination of features with functional groups
    if (len(core_features) >= 2) or (score >= 4 and len(all_features) >= 2):
        features_str = ", ".join(all_features)
        return True, f"Monoterpenoid features found: {features_str}"
    else:
        return False, "Insufficient characteristic features for monoterpenoid classification"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25409',
                          'name': 'monoterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'monoterpene. The term includes '
                                        'compounds in which the C10 skeleton '
                                        'of the parent monoterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['LIPID_MAPS_class:LMPR0102'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7499999999999999 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'C[C@]12[C@@H]3CC[C@@H](O3)[C@@]1(C)C(=O)OC2=O NAME: '
               'cantharidin REASON: CORRECT Monoterpenoid features found: '
               'cyclohexane, oxygen-containing groups, 3 ring(s)\n'
               ' * SMILES: C1CO[C@@](C[C@@H]1C)(C=C(C)C)[H] NAME: (2S,4R)-rose '
               'oxide REASON: CORRECT Monoterpenoid features found: gem '
               'dimethyl, isoprene unit, 1 ring(s)\n'
               ' * SMILES: CC(C)c1ccc(C)cc1O NAME: thymol REASON: CORRECT '
               'Monoterpenoid features found: gem dimethyl, oxygen-containing '
               'groups, 1 ring(s)\n'
               ' * SMILES: CC(C)C1(O)CCC(C)=CC1 NAME: 4-terpineol REASON: '
               'CORRECT Monoterpenoid features found: gem dimethyl, isoprene '
               'unit, oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: '
               'O=C1C=C(C)[C@]2([C@]1(O)[C@@]3(C(C)(C)C[C@@]2(O)C3)C)O NAME: '
               'Hitoyol A REASON: CORRECT Monoterpenoid features found: gem '
               'dimethyl, cyclohexane, cyclopentane, isoprene unit, '
               'oxygen-containing groups, 3 ring(s)\n'
               ' * SMILES: O[C@]12[C@]([C@H](C(=C)CC1)C2)(C(O)/C=C/C(O)(C)C)C '
               'NAME: E-beta-trans-5,8,11-trihydroxybergamot-9-ene REASON: '
               'CORRECT Monoterpenoid features found: gem dimethyl, '
               'cyclohexane, isoprene unit, oxygen-containing groups, 3 '
               'ring(s)\n'
               ' * SMILES: '
               'O=C1[C@H]2[C@H]3[C@@](C)(CC1)C[C@@H]2[C@@H](C(O)(C)C)C3 NAME: '
               'Norpupukeanane A REASON: CORRECT Monoterpenoid features found: '
               'gem dimethyl, cyclohexane, cyclopentane, menthane skeleton, '
               'oxygen-containing groups, 3 ring(s)\n'
               ' * SMILES: CC(O)\\C=C\\[C@]1(O)C(C)=CC(=O)CC1(C)C NAME: '
               '(6R)-vomifoliol REASON: CORRECT Monoterpenoid features found: '
               'gem dimethyl, isoprene unit, oxygen-containing groups, 1 '
               'ring(s)\n'
               ' * SMILES: CC(C)C12CC1C(C)C(=O)C2 NAME: thujone REASON: '
               'CORRECT Monoterpenoid features found: gem dimethyl, '
               'cyclohexane, cyclopentane, menthane skeleton, carane skeleton, '
               'oxygen-containing groups, 2 ring(s)\n'
               ' * SMILES: O=C1C(OC)=C(C(=CO)C([C@@]1(CC)C)=O)[C@H](CC)C NAME: '
               'Mycosphine C REASON: CORRECT Monoterpenoid features found: '
               'isoprene unit, oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: C1(C(CC(CC1C(C)=C)C(C)=C)C)C NAME: '
               '1,2-Dimethyl-3,5-bis(1-methylethenyl)cyclohexane REASON: '
               'CORRECT Monoterpenoid features found: cyclohexane, isoprene '
               'unit, 1 ring(s)\n'
               ' * SMILES: O=C1C[C@H]([C@@](CC1)(C)C=C)C(C)=C NAME: Geijerone '
               'REASON: CORRECT Monoterpenoid features found: cyclohexane, '
               'isoprene unit, oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: OC1=C(CCC(C1=O)C)C(C)C NAME: Gamma-Diosphenol '
               'REASON: CORRECT Monoterpenoid features found: gem dimethyl, '
               'isoprene unit, oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: O[C@H]1[C@@]2(C[C@@](C(C2=O)(C)C)(C1)[H])C NAME: '
               '6-endo-Hydroxyfenchone REASON: CORRECT Monoterpenoid features '
               'found: gem dimethyl, cyclohexane, cyclopentane, '
               'oxygen-containing groups, 2 ring(s)\n'
               ' * SMILES: OC(=O)C[C@@H](CC/C=C(\\C)/C(O)=O)C NAME: '
               '2,6R-Dimethyl-2E-octen-1,8-dioic acid REASON: CORRECT '
               'Monoterpenoid features found: isoprene unit, oxygen-containing '
               'groups\n'
               ' * SMILES: C1CC2CC1C3C2SC4=C(C3C5=CC=CO5)SC(=O)N4 NAME: '
               'LSM-18934 REASON: CORRECT Monoterpenoid features found: '
               'cyclohexane, cyclopentane, oxygen-containing groups, 5 '
               'ring(s)\n'
               ' * SMILES: O/N=C/C1=CC[C@@H](C(=C)C)CC1 NAME: Perillartine '
               'REASON: CORRECT Monoterpenoid features found: isoprene unit, '
               'oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: O=C([C@H](C1=C(C(O)=C(C)C(=C1C)O)C=O)C)C NAME: '
               'Redoxcitrinin REASON: CORRECT Monoterpenoid features found: '
               'oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: C[C@@H]1CC[C@@H]([C@@H](O)C1)C(C)(C)O NAME: '
               '1beta,3alpha,4beta-p-menthane-3,8-diol REASON: CORRECT '
               'Monoterpenoid features found: gem dimethyl, cyclohexane, '
               'menthane skeleton, oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: O[C@]1([C@H]([C@H]([C@](O)(CCCO)C)CC1)C)C NAME: '
               '11,12,15-trinorcycloneran-3,7,10-triol REASON: CORRECT '
               'Monoterpenoid features found: cyclopentane, oxygen-containing '
               'groups, 1 ring(s)\n'
               ' * SMILES: O(C1C(C(C)C)C(O)C=C(C1=O)CO)C(=O)CC(C)C NAME: '
               '[5-Hydroxy-3-(hydroxymethyl)-2-oxo-6-propan-2-ylcyclohex-3-en-1-yl] '
               '3-methylbutanoate REASON: CORRECT Monoterpenoid features '
               'found: isopropyl, gem dimethyl, isoprene unit, '
               'oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: O(C(CCC=C(C)C)(C)C=C)C(=O)C(C)C NAME: Linalyl '
               'isobutyrate REASON: CORRECT Monoterpenoid features found: gem '
               'dimethyl, isoprene unit, oxygen-containing groups\n'
               ' * SMILES: [H][C@@]12C[C@@]1(C[C@@H](O)C2=C)C(C)C NAME: '
               '(-)-cis-sabinol REASON: CORRECT Monoterpenoid features found: '
               'gem dimethyl, cyclohexane, cyclopentane, isoprene unit, carane '
               'skeleton, oxygen-containing groups, 2 ring(s)\n'
               ' * SMILES: O=CC1=C(O)C(=C(O)C=C1[C@H](CCCC)C)CC NAME: '
               'Anguillosporal REASON: CORRECT Monoterpenoid features found: '
               'oxygen-containing groups, 1 ring(s)\n'
               ' * SMILES: CC1(C)O[C@]2(C)CC[C@H]1C[C@@H]2O NAME: '
               '2-exo-hydroxy-1,8-cineole REASON: CORRECT Monoterpenoid '
               'features found: gem dimethyl, cyclohexane, menthane skeleton, '
               'oxygen-containing groups, 3 ring(s)\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O=C(OC)/C=C(/OC)\\C(\\OC)=C/C(=C\\C=C(\\[C@H](C[C@@]1(O[C@H]1C=C)C)C)/C)/C '
               'NAME: Haliangicin C REASON: MISSED Carbon count (22) outside '
               'typical range for monoterpenoids (7-15)\n'
               ' * SMILES: '
               'O=C1NCC2=C1C=C(O)C3=C2O[C@]([C@@H](O)C3)(CC/C=C(/CO)\\CCC=C(C)C)C '
               'NAME: Stachybotrin A REASON: MISSED Carbon count (23) outside '
               'typical range for monoterpenoids (7-15)\n'
               ' * SMILES: '
               'C1([C@@]2(C[C@@]3([C@](C1)(O[C@@](O3)([C@@H]2COC(C=4C=CC=CC4)=O)[H])C)O)[H])=O '
               'NAME: paeoniflorigenone REASON: MISSED Carbon count (17) '
               'outside typical range for monoterpenoids (7-15)\n'
               ' * SMILES: '
               'O=C1NCC2=C1C=C(O)C3=C2O[C@](C=C3)(CC/C=C(\\CCC=C(C)C)/C)C '
               'NAME: Chartarutine H REASON: MISSED Carbon count (23) outside '
               'typical range for monoterpenoids (7-15)\n'
               ' * SMILES: COC1=CC(=O)C(CC=C)=CC1(OC)C(\\C)=C\\c1ccc2OCOc2c1 '
               'NAME: Futoquinol REASON: MISSED Carbon count (21) outside '
               'typical range for monoterpenoids (7-15)\n'
               ' * SMILES: '
               'O(C(=O)CCCCCCCCCCCCCCCCC)CC=1C(=C(O)C(=C(OC)C1)C/C=C(/CC(=O)C=C(C)C)\\C)C=O '
               'NAME: Hericenone D REASON: MISSED Carbon count (37) outside '
               'typical range for monoterpenoids (7-15)\n'
               ' * SMILES: CC1=CC=CC(=C1)C=C2C(=O)N=C(S2)NN=C3CC4CCC3(C4(C)C)C '
               'NAME: '
               '5-[(3-methylphenyl)methylidene]-2-[2-(4,7,7-trimethyl-3-bicyclo[2.2.1]heptanylidene)hydrazinyl]-4-thiazolone '
               'REASON: MISSED Carbon count (21) outside typical range for '
               'monoterpenoids (7-15)\n'
               ' * SMILES: '
               'O[C@@]1([C@@]2([C@]([C@H](O)C1)(C(=CO[C@H]2O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)C(OC)=O)[H])[H])C '
               'NAME: Shanzhiside methyl ester REASON: MISSED Carbon count '
               '(17) outside typical range for monoterpenoids (7-15)\n'
               ' * SMILES: C1C2(C(CC(C(C)C)(C1)O2)OCC3=C(C=CC=C3)C)C NAME: '
               '1-methyl-2-[(2-methylbenzyl)oxy]-4-(propan-2-yl)-7-oxabicyclo[2.2.1]heptane '
               'REASON: MISSED Carbon count (18) outside typical range for '
               'monoterpenoids (7-15)\n'
               ' * SMILES: '
               'OC1(C2C(CC1O)C(=COC2OC3OC(C(O)C(O)C3O)COC(=O)/C=C/C4=CC=C(O)C=C4)C(O)=O)C '
               'NAME: Lippioside I REASON: MISSED Carbon count (25) outside '
               'typical range for monoterpenoids (7-15)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'S(C(S)CC)CCC',
                                     'name': 'xi-1-(Propylthio)-1-propanethiol',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'C[C@@H]([NH3+])C(=O)N[C@H](CO)C([O-])=O',
                                     'name': 'D-alanyl-D-serine zwitterion',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'CCN1C=CN(CC)C1=C2C=3C=CC(=C(C=4C=CC(=C(C5=CC=C(C(=C6C=CC2=N6)C7=[N+](CC)C=CN7CC)[N-]5)C8=[N+](CC)C=CN8CC)N4)C9=[N+](CC)C=CN9CC)N3.[Cl-].[Cl-].[Cl-].[Cl-].[Cl-].[Mn+3]',
                                     'name': 'Aeol 10150',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'OC(=O)\\C=C/C(=O)CC(=O)C(O)=O',
                                     'name': '3-maleylpyruvic acid',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'N=1C(C)=CN=CC1C=C',
                                     'name': '2-Ethenyl-6-methylpyrazine',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'OC(CCCCCCCCCCCCCC)C(N)CO',
                                     'name': 'C17-sphinganine',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'PA(22:0/20:2(11Z,14Z))',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'O[C@@H](CCCCCCCCCCC(O)=O)CCCC#C',
                                     'name': '12S-Hydroxy-16-heptadecynoic '
                                             'acid',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'},
                                 {   'smiles': 'S=P(OCC)(OCC)[O-]',
                                     'name': 'O,O-diethyl thiophosphate',
                                     'reason': 'Insufficient characteristic '
                                               'features for monoterpenoid '
                                               'classification'}],
    'sample_false_negatives': [   {   'smiles': 'OC(=O)C[C@@H](CC/C=C(\\C)/C(O)=O)C',
                                      'name': '2,6R-Dimethyl-2E-octen-1,8-dioic '
                                              'acid',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'CCC(C)(CCCC(C)C=O)OC',
                                      'name': '6-methoxy-2,6-dimethyloctanal',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'O([C@]12[C@]([C@@H](OC=C1)O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)(C(=C[C@H]2O)CO)[H])[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)CO',
                                      'name': 'Rehmannioside D',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'O=C(/C(=C/[C@@H]([C@H](O)[C@H](C[C@H](CC)C)C)C)/C)CC',
                                      'name': 'Pestalotiotone B',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'C[C@@H]1CCc2c(C)coc2C1',
                                      'name': '(+)-menthofuran',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'C(=C/C(/C)=C/C)(\\C=C(/C)\\C=C\\C)/C',
                                      'name': '3,5,7-Trimethyl-2E,4E,6E,8E-decatetraene',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'C(/C=C(/C=C/CC)\\C)(=C\\C(\\C)=C\\C)\\C',
                                      'name': '3,5,7-Trimethyl-2E,4E,6E,8E-undecatetraene',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'CC(C)(Cl)[C@H](Br)CC[C@@](Cl)(CBr)C(Cl)=C',
                                      'name': '(+)-halomon',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'CC1CCc2c(C)coc2C1',
                                      'name': 'menthofuran',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'},
                                  {   'smiles': 'O=C(O)/C(=C/C(CC(CC)C)C)/C',
                                      'name': 'Xylarinic acid B',
                                      'reason': 'Insufficient characteristic '
                                                'features for monoterpenoid '
                                                'classification'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Main Issue - False Negatives:\n'
                 'The program missed several valid monoterpenoids mainly '
                 'because of the carbon count restriction (7-15). Looking at '
                 'the false negatives, many are complex monoterpenoid '
                 'derivatives with additional functional groups or conjugated '
                 'systems, pushing their carbon count higher.\n'
                 '\n'
                 '2. Pattern Analysis:\n'
                 '- Most false negatives have the basic monoterpenoid skeleton '
                 'but with additional decorations (glycosylation, aromatic '
                 'rings, extended chains)\n'
                 '- The core C10 skeleton is present but modified in these '
                 'structures\n'
                 '- Some compounds have nitrogen-containing groups which '
                 "weren't considered in the original code\n"
                 '\n'
                 '3. Suggested Improvements:\n'
                 '\n'
                 'a) Modify carbon count restriction:\n'
                 '- Instead of using total carbon count, we should try to '
                 'identify the core C10 skeleton\n'
                 '- Could use a more flexible carbon range (7-25) to '
                 'accommodate derivatives\n'
                 '- Or remove the strict carbon count restriction entirely and '
                 'rely more on structural features\n'
                 '\n'
                 'b) Add new structural patterns:\n'
                 '- Include nitrogen-containing functional groups\n'
                 '- Add patterns for common modifications (glycosylation, '
                 'aromatic substitutions)\n'
                 '- Include more complex ring systems\n'
                 '\n'
                 'c) Improve scoring system:\n'
                 '- Make it more flexible for derivatives\n'
                 '- Give more weight to core monoterpenoid features\n'
                 '- Consider relative proportions rather than absolute counts\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 915,
    'num_false_positives': 100,
    'num_true_negatives': 16,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.9014778325123153,
    'recall': 0.9838709677419355,
    'f1': 0.9408740359897171,
    'accuracy': 0.890057361376673,
    'negative_predictive_value': 0.5161290322580645}