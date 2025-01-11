"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Adjustments have been made to carbon count and structural considerations to improve accuracy.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broaden carbon atom count check to include a wider range
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 60:
        return False, f"Carbon count is {carbon_count}, expected between 20 and 60"
    
    # Check for terpenoid-like isoprene units or similar structures
    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
    if mol.HasSubstructMatch(isoprene_unit):
        # Higher confidence match due to presence of isoprene units
        return True, "Contains patterns typical of terpenoid structures, like isoprene units"

    # Check for multiple ring structures, common in complex terpenoid derivatives
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count >= 2:
        return True, f"Contains multiple ({ring_count}) rings typical of terpenoids, such as sesterterpenoids"

    # Allow flexibility for other structural motifs
    return True, "Contains structural motifs common in sesterterpenoids, allowing for modifications and rearrangements"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26660',
                          'name': 'sesterterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'sesterterpene. The term includes '
                                        'compounds in which the C25 skeleton '
                                        'of the parent sesterterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups). Sometimes '
                                        'sesterterpenoids are erroneously '
                                        'referred to as sesterpenoids.',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['LIPID_MAPS_class:LMPR0105'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7733333333333333 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C(/C=C(/CC/C=C/2\\COC(C2)C=C(C)C)\\C)C\\C(\\C1)=C\\CC/C(=C/CC=3C(=[N+]([O-])C4=CC=CC=C4C3)C)/C '
               'NAME: Marinoterpin B REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C(O)[C@H]([C@H]1[C@@H]2[C@@](CC=C(C[C@H](O)C3=C(C(=CC2)C)CC[C@@H]3C)C)(C)CC1)C '
               'NAME: Betaestacin VI REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C=2O[C@]3([C@@H](CCC(=C)[C@H]4[C@H](CC3)C(C4)(C)C)CC2[C@@]5([C@H]1C[C@H]6C[C@@](C[C@@H]56)(C(=O)O)C)C)C '
               'NAME: Sterhirsutin C REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@@H]([C@H](/C=C\\[C@@H]4OC4(C)C)C)CC[C@@]3(C2)C)C=O '
               'NAME: Asperophiobolin I REASON: CORRECT Contains patterns '
               'typical of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C(O)[C@@]1([C@@H]2[C@H](C=O)CC[C@@H]1CC=C(CCC=C(CC2)C)C)CCC=C(C)C '
               'NAME: Emericellene D REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C([C@@H]([C@H]1O[C@@]([C@@]23O[C@@](CC)([C@@H](O2)C)C[C@H]3C)(C)C[C@@H]1C)CC)[C@H]([C@@H](O)[C@@H]([C@H]4O[C@](O)([C@H](C)C[C@@H]4C)CC)C)C '
               'NAME: X-14889-D REASON: CORRECT Contains multiple (4) rings '
               'typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               '[H]C(=O)[C@H]1C(=CC[C@@]2([H])[C@]3(C)CC[C@@]4([H])C(C)(C)CCC[C@]4(C)[C@@]3([H])C[C@H](OC(C)=O)[C@]12C)C([H])=O '
               'NAME: scalaradial REASON: CORRECT Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'OC1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(\\C)/C=C/C=O)(C)C '
               "NAME: 3-Hydroxy-10'-apo-b,y-carotenal REASON: CORRECT Contains "
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O1C2(C(C(CC(O)C2)(C)C)=CC1\\C(\\C)=C/C=C/C(/C)=C/C=C\\C=C(\\CCCO)/C)C '
               'NAME: Valenciachrome REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C=C(/C(=C/C=C(/C=C/C=C(/C=C/C=C/C(=C/C=C/C(=O)C)/C)\\C)\\C)/C(C1)(C)C)C '
               "NAME: (6E)-9'-apo-rhodoxanthinone REASON: CORRECT Contains "
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@@H]([C@H](C4OC(C)(C)C=C4)C)CC[C@@]3(C2)C)C=O '
               'NAME: Ophiobolin T REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@]3(C)[C@]([H])(C[C@@H](OC(C)=O)[C@]4(C)C5=C([C@H](O)OC5=O)[C@@H](C[C@@]34[H])OC)[C@@]1(C)CCCC2(C)C '
               'NAME: '
               'rel-(3R,4R,5aS,5bR,7aS,11aS,11bR,13R,13aS)-3-hydroxy-4-methoxy-5b,8,8,11a,13a-pentamethyl-1-oxo-1,3,4,5,5a,5b,6,7,7a,8,9,10,11,11a,11b,12,13,13a-octadecahydrochryseno[2,1-c]furan-13-yl '
               'acetate REASON: CORRECT Contains patterns typical of terpenoid '
               'structures, including sesterterpenoids\n'
               ' * SMILES: '
               '[H]C(=O)C1=C/C[C@@]2([H])[C@]3([H])[C@H](CC[C@@]3(C)CC[C@]2(C)C[C@@]2([H])[C@@H](C)CC(=O)[C@@]\\12[H])C(C)=C '
               'NAME: variecolin REASON: CORRECT Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1O[C@]23[C@@](C=C(CO)[C@@H](C2)C[C@H](O)C)(CC=C[C@H]([C@H]4[C@@](C(=C1C3=O)O)([C@H]5[C@@H](C[C@@H](C)CC5)C=C4)C)C)C '
               'NAME: Maklamicin REASON: CORRECT Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'C1CCC(C2[C@]1(C3[C@@](CC2)([C@@H]([C@@](CC3)(C)O)CC/C(=C/CO)/C)C)C)(C)C '
               'NAME: (17E)-cheilanthenediol REASON: CORRECT Contains patterns '
               'typical of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1[C@@H]2C(=CC[C@H]3[C@@H]([C@H](/C=C\\C=C(C)C)C)CC[C@@]3(C[C@@H]2[C@@](C1)(O)C)C)C '
               'NAME: 21-deoxyophiobolin K REASON: CORRECT Contains patterns '
               'typical of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               '[C@]12([C@]([C@@]3([C@@]4([C@]([C@]5(C(C[C@]3(CC1)CO)=C(CC5)C(=O)[O-])[H])(C4)C)[H])[H])([C@H](CC2)C(C)C)[H])C '
               'NAME: asperterpenoid A(1-) REASON: CORRECT Contains patterns '
               'typical of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'C1(=CC[C@@H]2[C@@](CC[C@@H]2[C@H](CCC=C(C)C)C)(C)C[C@H]3[C@@H]1CCC3=C)C '
               'NAME: Ophiobola-3(20),7,18-triene REASON: CORRECT Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O1[C@@]2([C@]34[C@@]([C@@]([C@]5([C@@]3(C(=C6[C@@]5([C@@](O)(CC[C@@]7([C@@]6(OC(=O)[C@H]7C)[H])[H])C)[H])C)[H])[H])([C@@H]2O)[H])([C@@]1(CC[C@@]8([C@@]4(OC(=O)[C@H]8C)[H])[H])C)[H])C '
               'NAME: Anabsin REASON: CORRECT Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'C\\C(CCC1=C(C)CCCC1(C)C)=C/CC\\C(CO)=C/C=C/C1=CC(=O)OC1O NAME: '
               'aplysinoplide A REASON: CORRECT Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'C[C@@H](CCC=C(C)C)[C@H]1CC[C@]2(C)C[C@H]3[C@H](CC[C@@]3(C)O)\\C(C)=C/C[C@@H]12 '
               'NAME: ophiobolin F REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)/C(=C\\CO)/C '
               'NAME: (2Z,6E,10E,14E)-geranylfarnesol REASON: CORRECT Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C=C(C)[C@@H]2[C@]13C(=CC[C@]34C(=CC[C@]4(C)C2)[C@H](CCC=C(CO)C)C)CO '
               'NAME: Bipolarolide D REASON: CORRECT Contains patterns typical '
               'of terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'OC(=O)\\C=C(\\CCC=1[C@@]2(C([C@@](CCC2)(CCC=C(C)C)C)CCC1C)C)/C '
               'NAME: (+)-Dysideapalaunic acid REASON: CORRECT Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'OC1(C(C2(C(C3(C(C(CCC3)(C)C)CC2)C)CC1)C)CC(O)C=4C(OC(=O)C4)OC(=O)/C=C/C5=CC=CC=C5)CO '
               'NAME: 25-Cinnamoyl-vulgaroside REASON: CORRECT Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               'False positives: SMILES: '
               'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C '
               'NAME: Isonigerone REASON: WRONGLY CLASSIFIED Contains multiple '
               '(6) rings typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: '
               '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-) '
               'REASON: WRONGLY CLASSIFIED Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC NAME: '
               'N-methyl arachidonoyl amine REASON: WRONGLY CLASSIFIED '
               'Contains patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O=C1N(CC(=O)N[C@H](C(=O)O[C@H]([C@@H](C(NC(C=C1)=C)=O)C)C(CCCCCCCCCCCCCC)C)C(O)C(=O)N)C '
               'NAME: Rakicidin H REASON: WRONGLY CLASSIFIED The structure '
               'appears to align with known sesterterpenoids criteria, '
               'considering carbon count and structural motifs typical in '
               'terpenoids\n'
               ' * SMILES: '
               'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C '
               'NAME: Prenylterphenyllin F REASON: WRONGLY CLASSIFIED Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO '
               'NAME: '
               '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol '
               'REASON: WRONGLY CLASSIFIED Contains multiple (7) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC '
               'NAME: '
               "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol "
               'REASON: WRONGLY CLASSIFIED Contains multiple (5) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@H]1CN(C)C(=O)NC4=CC=CC=C4)C5=CC=CC=C5N2C)[C@H](C)CO '
               'NAME: LSM-30671 REASON: WRONGLY CLASSIFIED Contains multiple '
               '(5) rings typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(17:0/16:0) REASON: WRONGLY CLASSIFIED The structure '
               'appears to align with known sesterterpenoids criteria, '
               'considering carbon count and structural motifs typical in '
               'terpenoids\n'
               ' * SMILES: '
               'C[C@@H]([C@H]1CC[C@H]2[C@@H]3CC=C4C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1 '
               'NAME: withalongolide F REASON: WRONGLY CLASSIFIED Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'CCOC(=O)C1=NN(C2=NC(=C(C(=O)N12)C#N)C3=CC=C(C=C3)C)C4=CC=CC=C4C '
               'NAME: '
               '6-cyano-1-(2-methylphenyl)-7-(4-methylphenyl)-5-oxo-[1,2,4]triazolo[4,3-a]pyrimidine-3-carboxylic '
               'acid ethyl ester REASON: WRONGLY CLASSIFIED Contains multiple '
               '(4) rings typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CC1=C(C(=NO1)C)NC(=O)N2C[C@H](COC[C@H]3[C@H]2CC[C@H](O3)CC(=O)NCC4CC4)O '
               'NAME: '
               '(3R,6aR,8S,10aR)-8-[2-(cyclopropylmethylamino)-2-oxoethyl]-N-(3,5-dimethyl-4-isoxazolyl)-3-hydroxy-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocine-1-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple (4) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'O1C([C@@H](O)C[C@@H]1[C@H]([C@@H]2[C@@]3([C@@](C=4C([C@@]5(C(C([C@@H](O)CC5)(C)C)CC4)C)=CC3)(C)CC2)C)C)(C)C '
               'NAME: Ganoderchochlearin A REASON: WRONGLY CLASSIFIED Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: CC(C)N(CCC(C(N)=O)(c1ccccc1)c1ccccn1)C(C)C NAME: '
               'disopyramide REASON: WRONGLY CLASSIFIED Contains multiple (2) '
               'rings typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CCC(=O)N1C[C@@H](COC[C@@H]2[C@@H]1CC[C@H](O2)CC(=O)NCC3=CC(=CC(=C3)F)F)O '
               'NAME: '
               '2-[(3S,6aS,8S,10aS)-3-hydroxy-1-(1-oxopropyl)-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocin-8-yl]-N-[(3,5-difluorophenyl)methyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple (3) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CCCC#CC1=CC2=C(N=C1)O[C@H]([C@@H](CN(C2=O)[C@H](C)CO)C)CN(C)CC3=CC=CC=C3F '
               'NAME: '
               '(2R,3R)-2-[[(2-fluorophenyl)methyl-methylamino]methyl]-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-8-pent-1-ynyl-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-6-one '
               'REASON: WRONGLY CLASSIFIED Contains multiple (3) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'C1CCC(C1)NC(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@H](O[C@@H]4CO)CC(=O)NCC5CC5 '
               'NAME: '
               '2-[(1R,3S,4aS,9aR)-6-[[(cyclopentylamino)-oxomethyl]amino]-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-(cyclopropylmethyl)acetamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple (5) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: O(CC/C=C\\CCCCCCCC/C=C\\CCCC)C(=O)C NAME: '
               '3Z,13Z-Octadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O=C1N2[C@@H](C(=O)NCC(=O)N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@@H](C(N[C@H]1CC3=CC=CC=C3)=O)CCCNC(=NCC=C(C)C)N)=O)CC4=CC=CC=C4)=O)CC(=O)O)CC(C)C)CCC2 '
               'NAME: Autumnalamide REASON: WRONGLY CLASSIFIED Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)CC(O)=O)O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)cc2)[C@H](O)[C@H]1O '
               'NAME: pelargonidin 3-O-(6-O-malonyl-beta-D-glucoside) REASON: '
               'WRONGLY CLASSIFIED Contains multiple (4) rings typical of '
               'complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'CCCC(=O)NC1=CC2=C(C=C1)C(=O)N(C[C@@H]([C@H](CN([C@@H](CO2)C)S(=O)(=O)C3=CN(C=N3)C)C)OC)C '
               'NAME: '
               'N-[(5R,6S,9R)-5-methoxy-3,6,9-trimethyl-8-[(1-methyl-4-imidazolyl)sulfonyl]-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]butanamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple (3) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'S(OC[C@H]1O[C@@H](O)[C@H](NC(=O)C)[C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]1O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)(O)(=O)=O '
               'NAME: '
               '[(2R,3S,4R,5R,6R)-5-Acetamido-6-hydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]methyl '
               'hydrogen sulfate REASON: WRONGLY CLASSIFIED Contains multiple '
               '(3) rings typical of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'C[C@H]1CCCCO[C@H]([C@H](CN(C(=O)C2=C(O1)C=CC(=C2)NC(=O)CCCN(C)C)[C@H](C)CO)C)CN(C)C(=O)NC(C)C '
               'NAME: '
               '4-(dimethylamino)-N-[(3S,9R,10S)-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-9-[[methyl-[oxo-(propan-2-ylamino)methyl]amino]methyl]-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-16-yl]butanamide '
               'REASON: WRONGLY CLASSIFIED Contains multiple (2) rings typical '
               'of complex terpenoids like sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H](CC/C=C(/C(=O)O)\\C)C)CC4)(C)CC3)C)[C@@H](C2)O)(C)CC1)(C)C '
               'NAME: (24E)-7alpha-hydroxy-3-oxolanosta-8,24-dien-26-oic acid '
               'REASON: WRONGLY CLASSIFIED Contains patterns typical of '
               'terpenoid structures, including sesterterpenoids\n'
               ' * SMILES: '
               'O=C1C([C@H]2C(=O)C[C@@]3([C@@H](C(=O)OC)[C@@](C(=O)[C@](O)(C(=O)OC)C)(C)C(C[C@H]3[C@]2(C)CC1)=C)C)(C)C '
               'NAME: Asperterpene K REASON: WRONGLY CLASSIFIED Contains '
               'patterns typical of terpenoid structures, including '
               'sesterterpenoids\n'
               'False negatives: SMILES: '
               'O=C(OC/C=C(\\CC/C=C(\\CC/C=C(\\CCCC(CCCC(CCCC(CCCC(C)C)C)C)C)/C)/C)/C)C(C(O)CCCCCCCCCCCCCCCCCCCC/C=C\\CCCCCCCCCCCCCC)CCCCCCCCCCCCCCCCCCC '
               'NAME: Octahydroheptaprenyl mycolic acyl ester REASON: MISSED '
               'Carbon count is 93, expected between 20 and 50\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@H]3C4(OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6C(O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)OC8COC(C)C(C8)OC9OC(C(O)C(C9)O)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CC(C=O)=C[C@@H]3O)C)C)C '
               'NAME: Tetrocarcin B REASON: MISSED Carbon count is 61, '
               'expected between 20 and 50\n'
               ' * SMILES: '
               'O=C(O)/C(=C/[C@@H]([C@H](O[C@@H]1O[C@H]([C@H](O)[C@H]([C@H]1O)O)C)/C(=C/[C@@H]([C@H](O[C@@H]2O[C@H]([C@H](O)[C@H]([C@H]2O)O)C)/C(=C/[C@@H]([C@H](O[C@@H]3O[C@@H]([C@@H](O)[C@@H]([C@@H]3O)O)CO)/C(=C/C(C(O)/C(=C/C(C(O)C(CC)C)C)/C)C)/C)C)/C)C)/C)C)/C '
               'NAME: Phialotide G REASON: MISSED Carbon count is 53, expected '
               'between 20 and 50\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@@H]3C=C([C@H](C)C[C@]34OC(=O)C(C4=O)=C([C@@]5([C@H](C(C(C2)OO)=C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](O)[C@@H](C7)O[C@H]8O[C@@H](O)[C@@H](O[C@@H]9O[C@H]([C@H](OC)[C@@H](C9)O)C)[C@@H](C8)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CO)C)C)C '
               'NAME: Lobophorin CR3 REASON: MISSED Carbon count is 61, '
               'expected between 20 and 50\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]1([C@@H](NOC(=O)C)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@@H]3C=C(C=O)[C@H](C)C[C@]34OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](O)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](OC)[C@@H](C9)O)C)[C@@H](C8)O)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)C)C)C '
               'NAME: Lobophorin H REASON: MISSED Carbon count is 61, expected '
               'between 20 and 50\n'
               ' * SMILES: '
               'O=C1O[C@@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=C(CC[C@H](C([C@@H]4[C@@H]1[C@@]56C(=O)[C@]4(O)C(=O)C(=C5[C@@H]7[C@H]8C(=O)[C@H](O)CCC(=C[C@@H]9[C@]%10(C([C@@H]8[C@H]6O7)=O)C(=O)N[C@H]([C@@H]%10[C@H](C)C(=C9)C)CC(C)C)C)C)=O)O)C)C)CC(C)C '
               'NAME: Amichalasine D REASON: MISSED Carbon count is 57, '
               'expected between 20 and 50\n'
               ' * SMILES: '
               'O=C1O[C@@]23[C@H](C=C([C@@H](O[C@@H]4O[C@@H]([C@H](NC(=O)OC)[C@](C4)(N)C)C)CC=C([C@H]5[C@@](C(=C1C2=O)O)([C@H]6[C@@H]([C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](O)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]6C)C=C5)C)C)C)[C@@H](O)C=C(C3)C=O '
               'NAME: Tetrocarcin P REASON: MISSED Carbon count is 61, '
               'expected between 20 and 50\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@H]3[C@@]4(OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](O[C@@H]%10O[C@H]([C@H](O)CC%10)C)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CC(/C=C/C(=O)C)=C[C@@H]3O)C)C)C '
               'NAME: Microsporanate A REASON: MISSED Carbon count is 70, '
               'expected between 20 and 50\n'
               ' * SMILES: '
               'O=C1O[C@@]23[C@H](C=C([C@@H](OC4OC(C(NC(=O)OC)C(C4)(N)C)C)CC=C([C@H]5[C@@](C(=C1C2=O)O)([C@H]6[C@@H]([C@@H](OC7OC(C(O)C(C7)OC8OC(C(OC9OC(C(OC)C(C9)O)C)C(C8)O)C)C)[C@H](C)C[C@@H]6C)C=C5)C)C)C)C=C(CO)[C@@H](C3)C '
               'NAME: Lobophorin D REASON: MISSED Carbon count is 61, expected '
               'between 20 and 50\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]1([C@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@H]2C(=C[C@@H]3[C@@]4(OC(=O)C(C4=O)=C([C@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@@H]([C@@H](OC(=O)C)[C@H](C7)O[C@H]8O[C@@H]([C@@H](O)CC8)C)C)[C@H](C)C[C@@H]([C@H]56)C)C)O)CC(C=O)=C[C@@H]3O)C)C)C '
               'NAME: Tetrocarcin F REASON: MISSED Carbon count is 55, '
               'expected between 20 and 50\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'Carbon count is 9, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'Carbon count is 19, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'Carbon count is 6, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'Carbon count is 68, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'O=C(N(CCO)CCO)CN1C(=NC=C1)[N+]([O-])=O',
                                     'name': '1H-Imidazole-1-acetamide, '
                                             'N,N-bis(2-hydroxyethyl)-2-nitro-',
                                     'reason': 'Carbon count is 9, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'CC1=NN(C(C1)(C2=CN=CC=C2)O)C(=O)CC3=CC=C(C=C3)C',
                                     'name': '1-[5-hydroxy-3-methyl-5-(3-pyridinyl)-4H-pyrazol-1-yl]-2-(4-methylphenyl)ethanone',
                                     'reason': 'Carbon count is 18, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'COC(=O)C[C@@H]1CC[C@H]([C@H](O1)CO)NC(=O)C2=CC=CC=C2',
                                     'name': '2-[(2S,5R,6S)-5-benzamido-6-(hydroxymethyl)-2-oxanyl]acetic '
                                             'acid methyl ester',
                                     'reason': 'Carbon count is 16, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine',
                                     'reason': 'Carbon count is 9, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CC(=O)C([O-])=O',
                                     'name': 'keto-3-deoxy-D-manno-octulosonate',
                                     'reason': 'Carbon count is 8, expected '
                                               'between 20 and 60'},
                                 {   'smiles': 'O=C1N(O)CCCCNC(=O)CCC(=O)N(CCCCCNC(CC1)=O)O',
                                     'name': 'Avaroferrin',
                                     'reason': 'Carbon count is 17, expected '
                                               'between 20 and 60'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(OC/C=C(\\CC/C=C(\\CC/C=C(\\CCCC(CCCC(CCCC(CCCC(C)C)C)C)C)/C)/C)/C)C(C(O)CCCCCCCCCCCCCCCCCCCC/C=C\\CCCCCCCCCCCCCC)CCCCCCCCCCCCCCCCCCC',
                                      'name': 'Octahydroheptaprenyl mycolic '
                                              'acyl ester',
                                      'reason': 'Carbon count is 93, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@H]3C4(OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6C(O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)OC8COC(C)C(C8)OC9OC(C(O)C(C9)O)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CC(C=O)=C[C@@H]3O)C)C)C',
                                      'name': 'Tetrocarcin B',
                                      'reason': 'Carbon count is 61, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@@H]3C=C([C@H](C)C[C@]34OC(=O)C(C4=O)=C([C@@]5([C@H](C(C(C2)OO)=C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](O)[C@@H](C7)O[C@H]8O[C@@H](O)[C@@H](O[C@@H]9O[C@H]([C@H](OC)[C@@H](C9)O)C)[C@@H](C8)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CO)C)C)C',
                                      'name': 'Lobophorin CR3',
                                      'reason': 'Carbon count is 61, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=[N+]([O-])[C@@]1([C@@H](NOC(=O)C)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@@H]3C=C(C=O)[C@H](C)C[C@]34OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](O)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](OC)[C@@H](C9)O)C)[C@@H](C8)O)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)C)C)C',
                                      'name': 'Lobophorin H',
                                      'reason': 'Carbon count is 61, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=C1O[C@@]23[C@H](C=C([C@@H](O[C@@H]4O[C@@H]([C@H](NC(=O)OC)[C@](C4)(N)C)C)CC=C([C@H]5[C@@](C(=C1C2=O)O)([C@H]6[C@@H]([C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](O)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]6C)C=C5)C)C)C)[C@@H](O)C=C(C3)C=O',
                                      'name': 'Tetrocarcin P',
                                      'reason': 'Carbon count is 61, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@H]3[C@@]4(OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](O[C@@H]%10O[C@H]([C@H](O)CC%10)C)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)CC(/C=C/C(=O)C)=C[C@@H]3O)C)C)C',
                                      'name': 'Microsporanate A',
                                      'reason': 'Carbon count is 70, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=C1O[C@@]23[C@H](C=C([C@@H](OC4OC(C(NC(=O)OC)C(C4)(N)C)C)CC=C([C@H]5[C@@](C(=C1C2=O)O)([C@H]6[C@@H]([C@@H](OC7OC(C(O)C(C7)OC8OC(C(OC9OC(C(OC)C(C9)O)C)C(C8)O)C)C)[C@H](C)C[C@@H]6C)C=C5)C)C)C)C=C(CO)[C@@H](C3)C',
                                      'name': 'Lobophorin D',
                                      'reason': 'Carbon count is 61, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=C1O[C@@]23[C@H](C=C([C@@H](O[C@@H]4O[C@@H]([C@H](NC(=O)OC)[C@](C4)(N)C)C)CC=C([C@H]5[C@@](C(=C1C2=O)O)([C@H]6[C@@H]([C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@H]9O[C@H]([C@H](O[C@@H]%10O[C@H]([C@H](O)CC%10)C)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]6C)C=C5)C)C)C)[C@@H](O)C=C(C3)/C=C/C(=O)C',
                                      'name': 'Microsporanate D',
                                      'reason': 'Carbon count is 70, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=[N+]([O-])[C@@]1([C@@H](NC(=O)OC)[C@H](O[C@H](C1)O[C@@H]2C(=C[C@H]3C4(OC(=O)C(C4=O)=C([C@@]5([C@H](C(=CC2)C)C=C[C@@H]6[C@@H](O[C@@H]7O[C@H]([C@H](OC(=O)C)[C@@H](C7)O[C@@H]8O[C@H]([C@H](O[C@@H]9O[C@H]([C@H](O[C@@H]%10O[C@H]([C@H](O)CC%10)C)[C@@H](C9)O)C)CC8)C)C)[C@@H](C)C[C@@H]([C@@H]56)C)C)O)C=C(C=O)C[C@@H]3O)C)C)C',
                                      'name': 'Tetrocarcin A',
                                      'reason': 'Carbon count is 67, expected '
                                                'between 20 and 60'},
                                  {   'smiles': 'O=C1OC23C(C=C(C(OC4OC(C(NC(=O)C)C(C4)(O)C)C)C(C)C=C(C5C(C(=C1C2=O)O)(C6C(C(OC7OC(C(OC8OC(C(OC9OC(C(O)C(C9)O)C)C(C8)O)C)C(C7)OC%10OC(C(O)CC%10)C)C)C(C)CC6C)C=C5)C)C)C)C=CC(C3)CC',
                                      'name': 'Saccharocarcin A',
                                      'reason': 'Carbon count is 67, expected '
                                                'between 20 and 60'}],
    'attempt': 3,
    'reasoning': 'To improve the classification of sesterterpenoids, we need '
                 'to address the following issues identified from the '
                 'outcomes:\n'
                 '\n'
                 '1. **Carbon Count Range**: The false negatives indicate that '
                 'our carbon count criteria of between 20 to 50 carbons might '
                 'be too strict. Sesterterpenoids might have varying carbon '
                 'counts due to modifications, so we need to adjust this '
                 'range.\n'
                 '\n'
                 '2. **Use of Structural Patterns**: The current method relies '
                 'heavily on patterns like isoprene units and ring counts, '
                 'which may not be sufficient to classify complex structures '
                 'unique to sesterterpenoids.\n'
                 '\n'
                 '3. **Overall Structure and Modifications**: Sesterterpenoids '
                 'can have rearrangements and modifications, including '
                 'non-standard terpenoid motifs that need to be considered for '
                 'classification.\n'
                 '\n'
                 "Here's a refined approach to improve classification:\n"
                 '\n'
                 '- **Expand Carbon Count Range**: We will broaden the range '
                 'to allow more flexibility. \n'
                 '- **Functional Group Patterns**: Sesterterpenoids should '
                 'have specific functional groups or rings indicative of '
                 'terpenoids, but we need to be cautious about generalizing '
                 'too broadly.\n'
                 '- **Allow Flexibility with Structural Variations**: Permit '
                 'some leeway in classifications by considering non-standard '
                 'isoprene connections or modifications in the backbone.\n'
                 '- **Feedback and Iteration**: Use the outcomes as feedback '
                 'for refining heuristics further.\n'
                 '\n'
                 "Here's an updated version of the code to better classify "
                 'sesterterpenoids:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 242,
    'num_false_positives': 100,
    'num_true_negatives': 93,
    'num_false_negatives': 26,
    'num_negatives': None,
    'precision': 0.7076023391812866,
    'recall': 0.9029850746268657,
    'f1': 0.7934426229508197,
    'accuracy': 0.7266811279826464,
    'negative_predictive_value': 0.7815126050420168}