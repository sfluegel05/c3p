"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid has a C20 skeleton derived from a diterpene, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Broaden carbon count range based on diterpenoids' variability
    if c_count < 15 or c_count > 50:
        return False, f"Uncommon carbon count ({c_count}) for diterpenoids"
    
    # Allow for diverse structures including cyclic, bicyclic, etc.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1 and c_count < 20:
        return False, "Structures like cyclic or polycyclic are typical for diterpenoids"

    # Check for double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if double_bond_count < 1:  # May check minimum double bond occurrence for diterpenoids
        return False, "Few double bonds; uncommon for diterpenoids"

    # Look for common functional groups
    possible_groups = [
        "[CX4][OX2H]",  # hydroxyl
        "[CX3](=O)",    # carbonyl
        "[OX2][CX4]",   # ether
        "C1OC1",        # epoxide
        "[CH3]",        # methyl
        "[OX2H]",       # hydroxyl group
    ]

    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in possible_groups):
        return False, "Missing typical functional groups like hydroxyls, carbonyls, ethers, epoxides, or methyl groups"

    # Ensure presence of at least one chiral center
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, "Lack of chiral centers; uncommon for diterpenoids"

    return True, "Molecule matches diterpenoid characteristics with flexible carbon count, diverse structures, and typical functional groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23849',
                          'name': 'diterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'diterpene. The term includes '
                                        'compounds in which the C20 skeleton '
                                        'of the parent diterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['LIPID_MAPS_class:LMPR0104'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.8983739837398375 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[C@@]12([C@@](C([C@H](CC1)OC(C[C@@H](N(C)C)C3=CC=CC=C3)=O)=C)([C@@H]([C@@]4(CC(C(=C([C@H]([C@@H]2O)OC(=O)C)C4(C)C)C)=O)O)O)[H])C '
               'NAME: taxine B REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               '[H][C@]12CC\\C(C)=C\\C[C@@]3(C)C(=O)C[C@H](C(C)=C)[C@]3([H])CC[C@]1(C)O2 '
               'NAME: '
               '(1R,3E,7S,8S,11S,12S)-7,8-Epoxy-14-oxo-3,18-dolabelladiene '
               'REASON: CORRECT Molecule matches diterpenoid characteristics '
               'with flexible carbon count, ring structures, and typical '
               'functional groups\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@H]3C[C@]1(C(=O)C3=C)[C@@]1(O)OC[C@]22CC[C@H](OC(C)=O)[C@@](C)(COC(C)=O)[C@@]2([H])[C@@H]1O '
               'NAME: maoecrystal J REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'O[C@H]1CC(C(=C(C1=O)C)C#C\\C(\\C)=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C#CC=2C(C[C@H](O)C(=O)C2C)(C)C)/C)/C)\\C)(C)C '
               "NAME: 7,8,7',8'-Tetradehydroastaxanthin REASON: CORRECT "
               'Molecule matches diterpenoid characteristics with flexible '
               'carbon count, ring structures, and typical functional groups\n'
               ' * SMILES: '
               'CCN1C[C@]2(COC(=O)c3ccccc3NC(C)=O)CC[C@H](OC)[C@@]34[C@@H]2C[C@@](O)(C13)[C@@]1(O)C[C@H](OC)[C@H]2C[C@]4(O)[C@]1(O)[C@H]2OC '
               'NAME: Finaconitine REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'CCN1C[C@@]2(CC[C@@H](C34[C@@H]1C(CC23)[C@]56C4C[C@@H](C(C5)C(=C)[C@H]6O)O)O)C '
               'NAME: LSM-6375 REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'O=C1C(OC)=C(OC)[C@H](OC(=O)C)[C@@H]([C@H]1C)CC=C(CCC=C(C[C@H](O)C=C(C(=O)O)C)C)C '
               'NAME: Antroquinonol U REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'OC([C@]1([C@]2([C@@](CC1)(C)C=C[C@@H](C[C@@H](O)C=C(C[C@@H]2O)C)C)[H])[H])(C)C '
               'NAME: (-)-2,7-Dolabelladiene-6beta,10alpha,18-triol REASON: '
               'CORRECT Molecule matches diterpenoid characteristics with '
               'flexible carbon count, ring structures, and typical functional '
               'groups\n'
               ' * SMILES: C12=C(CCC[C@@]1(CC[C@H](C2)C(=C)CCC=C(C)C)C)C NAME: '
               'Clavulatriene A REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'O[C@]1([C@H](/C=C/[C@@](O)([C@H](O)CO)C)[C@@]2([C@H](C(CCC2)(C)C)[C@H]([C@@H]1O)O)C)C '
               'NAME: Sterebin Q1 REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'C\\C1=C/CC\\C(C)=C\\CC(C)(C)[C@@H](O)CC\\C(C)=C\\CC1 NAME: '
               '(1S,4E,8E,12E)-2,2,5,9,13-pentamethylcyclopentadeca-4,8,12-trien-1-ol '
               'REASON: CORRECT Molecule matches diterpenoid characteristics '
               'with flexible carbon count, ring structures, and typical '
               'functional groups\n'
               ' * SMILES: '
               'O=C1C=C(C(C)C)[C@@H]2[C@]1(CC[C@@]3([C@@H]2C[C@@H]4[C@H](CO)[C@H]([C@]3(O4)O)O)C)C '
               'NAME: Cyafricanin H REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'O1CC2(C3(C(C4(C(CC3)C(C(O)CC4)(C)C)C)CCC2C(C(OC(=O)C)CC5OC(=O)C(=C5)C)C)C)CC1=O '
               'NAME: Hovenidulcigenin A REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'O[C@]1([C@@H]2C=C3CC[C@@H]4[C@@]([C@]3(C2)CC1)(CC[C@H]([C@]4(CO)C)O)C)CO '
               'NAME: Aphidicolin A22 REASON: CORRECT Molecule matches '
               'diterpenoid characteristics with flexible carbon count, ring '
               'structures, and typical functional groups\n'
               ' * SMILES: '
               'O=CC1=CC2(O)[C@@]3(CC[C@@]4(C([C@H]3C[C@H]1O2)=C(C(C)C)CC4)C)C '
               'NAME: Cyathin B3 REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'O=C1C(O)=C2[C@@](C=CC[C@]2(CO)C)(C)C3=C1C[C@](C=C)(C)CC3 NAME: '
               'Libertellenone G REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               '[C@]12(C[C@H]3[C@](C)(C(CC(=C(CC1)C)C2(C)C)O)CCC4OC[C@@]34O)[H] '
               'NAME: oxo-taxadiene-4alpha,9alpha-diol REASON: CORRECT '
               'Molecule matches diterpenoid characteristics with flexible '
               'carbon count, ring structures, and typical functional groups\n'
               ' * SMILES: '
               'O=C(O)/C=C(/C=C/[C@@H]1O[C@]2(O[C@@H]([C@@H](C)CC2)C/C=C(/C=C/[C@H](O)[C@H](/C=C/C(=O)O)C)\\C)CC[C@]1(OC(=O)CCC(=O)OC)CCCC)\\C '
               "NAME: Reveromycin A 4'-methyl ester REASON: CORRECT Molecule "
               'matches diterpenoid characteristics with flexible carbon '
               'count, ring structures, and typical functional groups\n'
               ' * SMILES: '
               'C1CCC([C@]2([C@]1([C@@]3([C@@H](CC2)C)O[C@@](CC3)(C)C=C)C)[H])(C)C '
               'NAME: 9,13(R)-epoxylabd-14-ene REASON: CORRECT Molecule '
               'matches diterpenoid characteristics with flexible carbon '
               'count, ring structures, and typical functional groups\n'
               ' * SMILES: '
               '[H][C@@]12C[C@@H](CC[C@@]1([H])[C@]1(C)[C@@H](O)CC[C@@](C)(COC(C)=O)[C@@]1([H])[C@H](O)C2=O)C(=C)C=O '
               'NAME: Eriocasin D REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'C[C@]12CCC[C@@]([C@@H]1CCC3=C[C@@](C(C[C@H]23)O)(C)C=C)(C)C(=O)O '
               'NAME: '
               '(1R,4aR,4bS,7R,10aR)-7-ethenyl-6-hydroxy-1,4a,7-trimethyl-3,4,4b,5,6,9,10,10a-octahydro-2H-phenanthrene-1-carboxylic '
               'acid REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'OC1C2C(C(C(=C1)C=O)CC/C(/COC(=O)C)=C\\C(O)=O)(CCCC2(C)C)C '
               'NAME: '
               '(E)-3-(Acetyloxymethyl)-5-(2-formyl-4-hydroxy-5,5,8a-trimethyl-1,4,4a,6,7,8-hexahydronaphthalen-1-yl)pent-2-enoic '
               'acid REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               '[H][C@]12CC[C@@]3([H])[C@@]45CO[C@](O)([C@@H](O)[C@]4([H])C(C)(C)CC[C@@H]5OC(C)=O)[C@]3([C@@H]1O)C(=O)C2=C '
               'NAME: Lasiokaurin REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: '
               'O=C1OC(=C(C)C2=C1C[C@@]3(C(=C[C@H](O)[C@H]4[C@]3(CC[C@@H]5[C@@]4(C=CC(=O)C5(C)C)C)C)C)O2)C '
               'NAME: Brevione I REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               ' * SMILES: OCC=1C2(C(C(C(CC2)C)(CCC(CO)=CC(O)=O)C)CCC1)C NAME: '
               '3-(Hydroxymethyl)-5-[5-(hydroxymethyl)-1,2,4a-trimethyl-2,3,4,7,8,8a-hexahydronaphthalen-1-yl]pent-2-enoic '
               'acid REASON: CORRECT Molecule matches diterpenoid '
               'characteristics with flexible carbon count, ring structures, '
               'and typical functional groups\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               '[H][C@]12[C@H](OC(=O)c3ccccc3)[C@]3(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c4ccccc4)c4ccccc4)C(C)=C([C@@H](OC(C)=O)C(=O)[C@]1(C)[C@@H](O)C[C@H]1OC[C@@]21OC(C)=O)C3(C)C '
               'NAME: paclitaxel REASON: MISSED Uncommon carbon count (47) for '
               'diterpenoids\n'
               ' * SMILES: OC1=C(C2=C(C3(C(C(CCC3)(C)C)CC2)C)C=C1)C(C)C NAME: '
               '(+)-Totarol REASON: MISSED Few double bonds; uncommon for '
               'diterpenoids\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](O)[C@@H](O[C@@H](OC/C(=C\\CC/C(=C/CC\\C(=C\\CC[C@@](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)CO)(C)C=C)\\C)/C)/C)[C@@H]1O)C)[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO[C@@H]5O[C@H]([C@H](O)[C@@H](O)[C@H]5O)C '
               'NAME: Capsianoside II REASON: MISSED Uncommon carbon count '
               '(50) for diterpenoids\n'
               ' * SMILES: '
               'O[C@@]12[C@@]3([C@H]([C@]4(CC[C@@H]3[C@H]4CC1)C)C[C@H]5[C@H](C(C)C)CC[C@@]5(C[C@@H]2C)C)C '
               'NAME: Sestermobaraol REASON: MISSED Few double bonds; uncommon '
               'for diterpenoids\n'
               ' * SMILES: '
               'CC(\\C=C\\C1=C(C)CCCC1(C)C)=C/C=C/C(C)=C/c1cc(\\C=C\\C=C(C)\\C=C\\C2=C(C)CCCC2(C)C)cc[n+]1CCO '
               'NAME: N-retinylidene-N-retinylethanolamine REASON: MISSED '
               'Uncommon carbon count (42) for diterpenoids\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(C5C(CC4)(CCC(C5)=C)C(OC6OC(C(O)C(O)C6O)CO)=O)=CC3)C)(CC2)C)(CC1)C)(C)C)C7OC(C(O)C(OC8OCC(O)C(O)C8O)C7O)C(O)=O '
               'NAME: '
               '6-{[4,4,6a,6b,14b-pentamethyl-11-methylidene-8a-({[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy}carbonyl)-1,2,3,4,4a,5,6,6a,6b,7,8,8a,9,10,11,12,12a,14,14a,14b-icosahydropicen-3-yl]oxy}-3,5-dihydroxy-4-[(3,4,5-trihydroxyoxan-2-yl)oxy]oxane-2-carboxylic '
               'acid REASON: MISSED Uncommon carbon count (46) for '
               'diterpenoids\n'
               ' * SMILES: '
               'CO[C@@H]1OC[C@@]2(C)CCC[C@]11[C@H]2CCc2cc(C(C)C)c(O)cc12 NAME: '
               '(-)-(4S,5S,10R,20R)-12,18-dihydroxyabieta-8,11,13-trien-20-aldehyde18,20-methyl '
               'acetal REASON: MISSED Few double bonds; uncommon for '
               'diterpenoids\n'
               ' * SMILES: '
               'O=C1[C@@H]([C@@]23[C@@H](C(=O)O)CC[C@@H]([C@@H]2C1)C(C3)(C)C)CO '
               'NAME: (+)-5(6)-dihydro-6-hydroxyterrecyclic acid A REASON: '
               'MISSED Uncommon carbon count (15) for diterpenoids\n'
               ' * SMILES: CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO NAME: '
               '(E,E,E)-geranylgeraniol REASON: MISSED Diterpenoids typically '
               'have ring structures\n'
               ' * SMILES: '
               'O1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@]3([C@]4(C([C@@]5(C(CC4)C([C@@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO[C@@H]7O[C@H]([C@H](O)[C@@H](O)[C@H]7O)C)CC5)(C)C)C)CC[C@@]3([C@]8([C@@]1([C@H](OC8=O)CC=C(C)C)C)[H])[H])C)C '
               'NAME: Mabioside A REASON: MISSED Uncommon carbon count (48) '
               'for diterpenoids\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'OC1=C(C2=C(C3(C(C(CCC3)(C)C)CC2)C)C=C1)C(C)C',
                                      'name': '(+)-Totarol',
                                      'reason': 'Few double bonds; uncommon '
                                                'for diterpenoids'},
                                  {   'smiles': 'O[C@@]12[C@@]3([C@H]([C@]4(CC[C@@H]3[C@H]4CC1)C)C[C@H]5[C@H](C(C)C)CC[C@@]5(C[C@@H]2C)C)C',
                                      'name': 'Sestermobaraol',
                                      'reason': 'Few double bonds; uncommon '
                                                'for diterpenoids'},
                                  {   'smiles': 'CC(\\C=C\\C1=C(C)CCCC1(C)C)=C/C=C/C(C)=C/c1cc(\\C=C\\C=C(C)\\C=C\\C2=C(C)CCCC2(C)C)cc[n+]1CCO',
                                      'name': 'N-retinylidene-N-retinylethanolamine',
                                      'reason': 'Lack of chiral centers; '
                                                'uncommon for diterpenoids'},
                                  {   'smiles': 'CO[C@@H]1OC[C@@]2(C)CCC[C@]11[C@H]2CCc2cc(C(C)C)c(O)cc12',
                                      'name': '(-)-(4S,5S,10R,20R)-12,18-dihydroxyabieta-8,11,13-trien-20-aldehyde18,20-methyl '
                                              'acetal',
                                      'reason': 'Few double bonds; uncommon '
                                                'for diterpenoids'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO',
                                      'name': '(E,E,E)-geranylgeraniol',
                                      'reason': 'Lack of chiral centers; '
                                                'uncommon for diterpenoids'},
                                  {   'smiles': 'O=C(OC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C=C\\C(=C/C=C/C(=C\\C=C\\C=C(\\C=C/C=C(\\C=C/C(=O)OC)/C)/C)/C)\\C',
                                      'name': 'trans-Geranylgeranylbixin',
                                      'reason': 'Lack of chiral centers; '
                                                'uncommon for diterpenoids'},
                                  {   'smiles': 'O=C(O)C1=C(C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C(=C(C(=O)O)C=C2)C)C)\\C)\\C)/C)/C)C=C1)C)C',
                                      'name': 'Synechoxanthin',
                                      'reason': 'Lack of chiral centers; '
                                                'uncommon for diterpenoids'},
                                  {   'smiles': 'C=12[C@]3(C)[C@@]([H])(CCC1C=C(C(=C2O)O)C(C)C)C(CCC3)(C)C',
                                      'name': '11-hydroxyferruginol',
                                      'reason': 'Few double bonds; uncommon '
                                                'for diterpenoids'},
                                  {   'smiles': 'OC/C=C(/C=C\\C=C(\\C=C\\C=1C(CC=CC1C)(C)C)/C)\\C',
                                      'name': '11-cis-3,4-didehydroretinol',
                                      'reason': 'Lack of chiral centers; '
                                                'uncommon for diterpenoids'},
                                  {   'smiles': 'CCN1C[C@]2(COC)CC[C@H](O)[C@@]34[C@@H]5C[C@H]6[C@H](OC)[C@@H]5[C@](O)(C[C@@H]6OC)[C@@](O)([C@@H](OC)[C@H]23)C14',
                                      'name': 'Delsoline',
                                      'reason': 'Few double bonds; uncommon '
                                                'for diterpenoids'}],
    'attempt': 5,
    'reasoning': "To improve our classification of diterpenoids, let's analyze "
                 'what went wrong in the previous attempt:\n'
                 '\n'
                 '### Observed Issues\n'
                 '1. **Carbon Count Check**: The constraint of having a carbon '
                 'count within 16 to 40 was too restrictive, causing some '
                 'diterpenoids with uncommon carbon skeletons to be missed.\n'
                 '2. **Ring Presence Check**: The requirement for ring '
                 'structures, though typical, was inflexible. There are '
                 'diterpenoids with linear or branch-like structures.\n'
                 '3. **Functional Groups**: Some diterpenoids might lack the '
                 'characteristic functional groups considered in the previous '
                 'function.\n'
                 '4. **Stereochemistry Check**: While many diterpenoids are '
                 'chiral, there might be metabolites with few or no chiral '
                 'centers.\n'
                 '\n'
                 '### Improvements\n'
                 '- **Adapt Carbon Count**: Broaden the acceptable range of '
                 'carbon count.\n'
                 '- **Adapt Ring Check**: Allow for both cylic and acyclic '
                 'structures, enabling diterpenoids with fewer or more complex '
                 'rings to be included.\n'
                 '- **Expand Pattern Library**: Add more functional groups '
                 'related to diterpenoids and ensure that at least some '
                 'stereo-centers are present.\n'
                 '\n'
                 '### Revised Code',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 752,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.8826291079812206,
    'f1': 0.9376558603491272,
    'accuracy': 0.8826291079812206,
    'negative_predictive_value': 0.0}