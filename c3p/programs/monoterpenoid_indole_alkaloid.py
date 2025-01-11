"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    These compounds are derived from tryptophan and typically secologanin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical MIA range"

    # Count basic statistics
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    ring_count = rdMolDescriptors.CalcNumRings(mol)

    # Basic requirements
    if n_count == 0:
        return False, "No nitrogen atoms found - required for alkaloid"
    if c_count < 15:
        return False, f"Too few carbons ({c_count}) for MIA structure"
    if ring_count < 3:
        return False, f"Too few rings ({ring_count}) for MIA structure"

    # Look for indole or modified indole cores
    indole_patterns = [
        "c1ccc2[nH]ccc2c1",  # Basic indole
        "c1ccc2nccc2c1",     # Modified indole
        "c1ccc2N=CCc2c1",    # Another modified form
        "c1ccc2NCCc2c1",     # Dihydroindole
        "C1=CC=C2C(=C1)NC=C2",  # Alternative representation
        "C1=CC=C2C(=C1)N=CC2",  # Another variant
        "[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]:1:[#7]:[#6]:[#6]:2"  # Generic form
    ]
    
    has_indole = False
    for pattern in indole_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_indole = True
            break
            
    if not has_indole:
        return False, "No indole or modified indole core found"

    # Look for characteristic MIA structural features
    mia_scaffold_patterns = [
        # Common MIA ring fusion patterns
        "[#6]1[#6]2[#7][#6][#6]1[#6][#6]2",  # Basic aspidosperma type
        "[#6]1[#6]2[#7][#6][#6]1[#6][#6]2[#8,#7]",  # Oxygenated variant
        "[#6]1[#6]2[#7][#6][#6]([#6]1)[#6][#6]2",  # Strychnos type
        "[#6]1[#6]2[#7][#6][#6]3[#6]1[#6][#6]2[#6]3",  # Iboga type
        "[#6]1[#6]2[#7][#6][#6]([#6]1)[#6][#6]2[#6](=O)",  # Carbonyl variant
        # Bridge patterns
        "[#6]1[#6]2[#6][#7][#6]1[#6][#6]2",  # Common bridge
        "[#6]1[#6]2[#6][#7][#6]1[#6][#6]2[#8]",  # Oxygenated bridge
    ]
    
    mia_pattern_matches = 0
    for pattern in mia_scaffold_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            mia_pattern_matches += 1

    # Check for secologanin-derived features
    secologanin_patterns = [
        "[#6]-[#6](=O)-O[#6]",  # Ester group
        "[#6]=[#6]-[#6]",  # Vinyl group
        "[#6]-[#6](O)-[#6]",  # Alcohol
        "[#6]1-[#6]-[#6]-[#8]-[#6]-1",  # Pyran ring
    ]
    
    secologanin_matches = 0
    for pattern in secologanin_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            secologanin_matches += 1

    # Calculate revised complexity score
    complexity_contributors = [
        ring_count * 2,          # Rings
        n_count * 2,             # Nitrogens
        mia_pattern_matches * 5,  # MIA patterns (increased weight)
        secologanin_matches * 3,  # Secologanin features
        (1 if o_count > 0 else 0) * 2  # Presence of oxygen
    ]
    complexity_score = sum(complexity_contributors)
    
    if complexity_score < 12:
        return False, f"Complexity score ({complexity_score}) too low for typical MIA"

    # Additional structural requirements
    if ring_count > 12 and mia_pattern_matches == 0:
        return False, "Complex structure lacks characteristic MIA patterns"

    # Check C/N ratio (typically between 6-25 for MIAs)
    cn_ratio = c_count / n_count
    if cn_ratio < 6 or cn_ratio > 25:
        return False, f"C/N ratio ({cn_ratio:.1f}) outside typical MIA range"

    features = []
    if mia_pattern_matches > 0:
        features.append(f"{mia_pattern_matches} MIA scaffold patterns")
    if secologanin_matches > 0:
        features.append(f"{secologanin_matches} secologanin-derived features")
    
    feature_str = ", ".join(features) if features else "characteristic structural features"
    return True, f"Contains indole core, complex ring system ({ring_count} rings), and {feature_str}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65323',
                          'name': 'monoterpenoid indole alkaloid',
                          'definition': 'A terpenoid indole alkaloid which is '
                                        'biosynthesised from L-tryptophan and '
                                        'diisoprenoid (usually secolaganin) '
                                        'building blocks.',
                          'parents': ['CHEBI:65321'],
                          'xrefs': [   'PMID:18280746',
                                       'PMID:20717879',
                                       'PMID:21425787',
                                       'PMID:22679912'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.4922279792746114 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[H][C@@]12N3CC[C@]11C(Nc4ccccc14)=C(C[C@]2(CC)C=CC3)C(=O)OC '
               'NAME: tabersonine REASON: CORRECT Contains indole core, '
               'complex ring system (5 rings), and ester group, ethyl/vinyl '
               'group, 3 sp2 carbons\n'
               ' * SMILES: '
               'C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC=CC6)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O '
               'NAME: pseudostrychnine REASON: CORRECT Contains indole core, '
               'complex ring system (7 rings), and ethyl/vinyl group, 2 sp2 '
               'carbons\n'
               ' * SMILES: '
               'N1[C@@]2([C@@]3(C4=C1C=CC=C4)[C@@H](C5[C@@H]([N+]6([C@]2(C[C@H]5[C@@H]([C@H]6O)CC)[H])C)C3)O)[H] '
               'NAME: 4-methylnorajmaline REASON: CORRECT Contains indole '
               'core, complex ring system (7 rings), and ethyl/vinyl group\n'
               ' * SMILES: '
               'CCC1(CC2CC(C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)C78CCN9C7C(C=CC9)(C(C(C8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O '
               'NAME: LSM-36946 REASON: CORRECT Contains indole core, complex '
               'ring system (9 rings), and ester group, ethyl/vinyl group, 4 '
               'sp2 carbons\n'
               ' * SMILES: '
               'N1([C@@]2([C@@]3(C4=C1C=CC=C4)[C@@H](C5[C@@H]([N+]6([C@]2(C[C@H]5[C@@H]([C@H]6O)CC)[H])C)C3)O)[H])C '
               'NAME: 4-methylajmaline REASON: CORRECT Contains indole core, '
               'complex ring system (7 rings), and ethyl/vinyl group\n'
               ' * SMILES: '
               'CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@](C=CC9)([C@H]([C@@]([C@@H]8N6C=O)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O '
               'NAME: LSM-15607 REASON: CORRECT Contains indole core, complex '
               'ring system (9 rings), and ester group, ethyl/vinyl group, 5 '
               'sp2 carbons\n'
               ' * SMILES: '
               'C/C=C\\1/C[N+]2(CC[C@@]34CC=5C=C6C(=CC5O[C@@]74[C@@]2(C[C@@]1([C@@](C(=O)OC)(N7C8=CC=CC=C83)[H])[H])[H])N[C@]9%10CC[C@@]%11(CC)CCC[N@+]%10(CC[C@@]69O)C%11)[O-] '
               'NAME: goniomedine A N-oxide REASON: CORRECT Contains indole '
               'core, complex ring system (11 rings), and ester group, '
               'ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'CC[C@]1(CN2CC[C@@]34C5=CC=CC=C5N6C[C@@]7(C[C@@]8(C[C@](CC)(N9CC[C@@]%10(C%11=CC=CC=C%11N=C7%10)C89)[H])[H])CC([C@]1(CC32)[H])=C46)[H] '
               'NAME: Leucoridine A REASON: CORRECT Contains indole core, '
               'complex ring system (11 rings), and ethyl/vinyl group, 2 sp2 '
               'carbons\n'
               ' * SMILES: '
               'C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC(=C(C6)OC)OC)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O '
               'NAME: pseudobrucine REASON: CORRECT Contains indole core, '
               'complex ring system (7 rings), and ethyl/vinyl group, 2 sp2 '
               'carbons\n'
               ' * SMILES: '
               'CC[C@@]12CC[C@@]34[C@@](CCN(C[C@@]5([C@]1(O5)[H])[H])[C@@]23[H])(C=6C=C(C=CC6N4)[C@]7([C@@]8([C@]([C@@]9(CC)CC(=C%10[C@]%11(CCN7[C@@]9%11[H])C%12=CC=CC=C%12N%10)C(=O)OC)(O8)[H])[H])[H])[H] '
               'NAME: Criophylline REASON: CORRECT Contains indole core, '
               'complex ring system (12 rings), and ester group, ethyl/vinyl '
               'group, 2 sp2 carbons\n'
               ' * SMILES: '
               'C/C=C\\1/CN2CC[C@@]34CC=5C=C6C=7CCN8CCC[C@](CC)(CCC7NC6=CC5O[C@@]94[C@]2(C[C@]1([C@@](C(=O)OC)(N9C%10=CC=CC=C%103)[H])[H])[H])C8 '
               'NAME: Goniomedine B REASON: CORRECT Contains indole core, '
               'complex ring system (10 rings), and ester group, ethyl/vinyl '
               'group, 2 sp2 carbons\n'
               ' * SMILES: '
               'C/C=C\\1/CN2CC[C@@]3(C4=CC=CC=C4N5C=C([C@]1(C[C@]2([C@]35O)[H])[H])C(=O)OC)O '
               'NAME: 2,7-dihydroxyapogeissoschizine REASON: CORRECT Contains '
               'indole core, complex ring system (5 rings), and ester group, '
               'ethyl/vinyl group, 3 sp2 carbons\n'
               ' * SMILES: '
               'C=12[C@@]3([C@H]4N(CC3)C\\C(\\[C@@](C1C(=O)OC)(C4)[H])=C\\C)C5=CC=CC=C5N2 '
               'NAME: akuammicine REASON: CORRECT Contains indole core, '
               'complex ring system (5 rings), and ester group, ethyl/vinyl '
               'group, 3 sp2 carbons\n'
               ' * SMILES: '
               '[H][C@@]12Nc3cc(OC)ccc3[C@@]11CCN3CC=C[C@@](CC)(C[C@]2(O)C(=O)OC)[C@@]13[H] '
               'NAME: 3-hydroxy-16-methoxy-2,3-dihydrotabersonine REASON: '
               'CORRECT Contains indole core, complex ring system (5 rings), '
               'and ester group, ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'C12=C(C=CC=C1)N=C3[C@]([C@@]4(/C(/C[NH+]5[C@@H](C4)[C@]32CC5)=C\\C)[H])(C(=O)OC)C(O)O '
               'NAME: 17-dehydropreakuammicine hydrate(1+) REASON: CORRECT '
               'Contains indole core, complex ring system (5 rings), and ester '
               'group, ethyl/vinyl group, 3 sp2 carbons\n'
               ' * SMILES: '
               'COC([C@@H]1C[C@@]23CCC[N+]4(CCC5(C6=CC=CC=C6N(C(=O)OC)[C@]15CC2)[C@]34[H])[O-])=O '
               'NAME: Pleiocarpoline REASON: CORRECT Contains indole core, '
               'complex ring system (7 rings), and ester group, ethyl/vinyl '
               'group, 2 sp2 carbons\n'
               ' * SMILES: '
               'CC[C@@]1(O)C[C@H]2C[N@@](C1)CCc1c([nH]c3ccccc13)[C@@](C2)(C(=O)OC)c1cc2c(cc1OC)N(C)[C@@H]1[C@]22CCN3CC=C[C@](CC)([C@@H]23)[C@@H](OC(C)=O)[C@]1(O)C(=O)OC '
               'NAME: Leurosidine REASON: CORRECT Contains indole core, '
               'complex ring system (9 rings), and ester group, ethyl/vinyl '
               'group, 4 sp2 carbons\n'
               ' * SMILES: '
               'O=C1N2[C@@]3([C@@]4([C@]5(N(CC4)C\\C(\\[C@](C5)(C3=CC1)[H])=C\\CO)[H])C=6C2=CC=CC6)[H] '
               'NAME: isostrychnine REASON: CORRECT Contains indole core, '
               'complex ring system (6 rings), and ethyl/vinyl group, 3 sp2 '
               'carbons\n'
               ' * SMILES: '
               '[C@@]123[C@@H](NC4=C1C=CC=C4)[C@](C[C@]5([C@@H]2N(CC=C5)CC3)CC)(C(=O)OC)O '
               'NAME: (3R)-3-hydroxy-2,3-dihydrotabersonine REASON: CORRECT '
               'Contains indole core, complex ring system (5 rings), and ester '
               'group, ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'C/C=C\\1/CN2CC[C@]34CC5C6CC7(C8=C(CC(C6CO[C@]5(C)O[C@]94[C@@]2(C[C@@]1([C@@](C(=O)OC)(N9C%10=CC=CC=C%103)[H])[H])[H])N7C)C%11=CC=CC=C%11N8C)[H] '
               'NAME: Villalstonine REASON: CORRECT Contains indole core, '
               'complex ring system (11 rings), and ester group, ethyl/vinyl '
               'group, 2 sp2 carbons\n'
               ' * SMILES: '
               '[C@@]123[C@@](N(C4=C1C=CC(=C4)O)C)([C@]([C@@H]([C@]5([C@@]2(N(CC=C5)CC3)[H])CC)O)(C(=O)OC)O)[H] '
               'NAME: 11-O-demethyl-17-O-deacetylvindoline REASON: CORRECT '
               'Contains indole core, complex ring system (5 rings), and ester '
               'group, ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'CC(C(C(N1C2=CC=CC=C2[C@]34CCN5C[C@]6([C@@H](C)OC[C@]([C@]6(C[C@@]35O)[H])([C@@]41[H])[H])[H])=O)O)=O '
               'NAME: Acetyl-splendoline REASON: CORRECT Contains indole core, '
               'complex ring system (6 rings), and ethyl/vinyl group, 2 sp2 '
               'carbons\n'
               ' * SMILES: '
               'O=C1N2[C@@]3([C@@]4([C@]5([N+](CC4)(C\\C(\\[C@](C5)(C3=CC1)[H])=C\\CO)[O-])[H])C=6C2=CC=CC6)[H] '
               'NAME: isostrychnine N-oxide REASON: CORRECT Contains indole '
               'core, complex ring system (6 rings), and ethyl/vinyl group, 3 '
               'sp2 carbons\n'
               ' * SMILES: '
               'CC[C@]12C[N@@]3C[C@@H](C[C@@](C(=O)OC)(c4[nH]c5ccccc5c4CC3)c3cc4c(cc3OC)N(C)[C@@H]3[C@]44CCN5CC=C[C@](CC)([C@@H]45)[C@@H](OC(C)=O)[C@]3(O)C(=O)OC)[C@H]1O2 '
               'NAME: Leurosine REASON: CORRECT Contains indole core, complex '
               'ring system (10 rings), and ester group, ethyl/vinyl group, 4 '
               'sp2 carbons\n'
               ' * SMILES: '
               'C/C=C\\1/CN2[C@]3(C[C@@]1(C4([C@@]2(CC5(C6=CC(=CC=C6N(C)[C@]35[H])O)[C@@]4(OC(\\C=C\\C7=CC(=C(C(=C7)OC)OC)OC)=O)[H])[H])C(=O)OC)[H])[H] '
               'NAME: O-trimethoxy-3,4,5-cinnamate-OH-vincamajine REASON: '
               'CORRECT Contains indole core, complex ring system (8 rings), '
               'and ester group, ethyl/vinyl group, 4 sp2 carbons\n'
               'False positives: SMILES: '
               '[H][C@@]12CCN3CC[C@@]4(C(Nc5ccccc45)=C1C(=O)OC)[C@@]3([H])[C@H]2CC '
               'NAME: 20(S)-tubotaiwine REASON: WRONGLY CLASSIFIED Contains '
               'indole core, complex ring system (5 rings), and ester group, '
               'ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'CCS(=O)(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)CC5=C(C=CC(=C5)F)F '
               'NAME: '
               "[(1R)-2-[(2,5-difluorophenyl)methyl]-1'-ethylsulfonyl-7-methoxy-9-methyl-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (5 rings), and ethyl/vinyl group\n'
               ' * SMILES: '
               'COC1=CC2=C(C=C1)C3=C(N2)[C@H](N(CC34CN(C4)C(=O)C5CCC5)CC6=C(C=CC(=C6)F)F)CO '
               'NAME: '
               "cyclobutyl-[(1S)-2-[(2,5-difluorophenyl)methyl]-1-(hydroxymethyl)-7-methoxy-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]yl]methanone "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (6 rings), and ethyl/vinyl group, 1 sp2 carbons\n'
               ' * SMILES: '
               '[H][C@@]12C[C@]3([H])N([C@H]4C[C@@]5([C@H](OC(C)=O)C14)c1ccccc1N[C@@]35[H])[C@H](O)\\C2=C\\C '
               'NAME: 1,2-dihydrovomilenine REASON: WRONGLY CLASSIFIED '
               'Contains indole core, complex ring system (7 rings), and ester '
               'group, ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'O=C/1NC(C(=O)N\\C1=C/C=2C3=C(C=CC(=C3)CC4OC4(C)C)NC2C(C=C)(C)C)C '
               'NAME: Epoxyisoechinulin A REASON: WRONGLY CLASSIFIED Contains '
               'indole core, complex ring system (4 rings), and ethyl/vinyl '
               'group, 4 sp2 carbons\n'
               ' * SMILES: '
               'OCC1[C@@]2(N(C=3C([C@@]42[C@]5(N(CC4)C\\C(\\[C@@]1(C5)[H])=C\\C)[H])=CC=CC3)C(=O)C)[H] '
               'NAME: '
               '1-[(1R,9S,11R,12E,17S)-12-Ethylidene-10-(hydroxymethyl)-8,14-diazapentacyclo[9.5.2.01,9.02,7.014,17]octadeca-2,4,6-trien-8-yl]ethanone '
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (5 rings), and ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'COC(=O)C1=C2Nc3ccccc3[C@@]22CCN3CC=C[C@](C1)(C(C)O)[C@@H]23 '
               'NAME: 19-Hydroxytabersonine REASON: WRONGLY CLASSIFIED '
               'Contains indole core, complex ring system (5 rings), and ester '
               'group, ethyl/vinyl group, 3 sp2 carbons\n'
               ' * SMILES: '
               'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)C=C4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C '
               'NAME: Teraspiridole A_130090 REASON: WRONGLY CLASSIFIED '
               'Contains indole core, complex ring system (7 rings), and ester '
               'group, ethyl/vinyl group, 4 sp2 carbons\n'
               ' * SMILES: '
               'CCS(=O)(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)CC5=CC6=C(C=C5)OCO6 '
               'NAME: '
               "[(1R)-2-(1,3-benzodioxol-5-ylmethyl)-1'-ethylsulfonyl-7-methoxy-9-methyl-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (6 rings), and ethyl/vinyl group\n'
               ' * SMILES: '
               '[H][C@]12C[C@]34C(OC[C@@]1([H])[C@]1([H])C[C@]3([H])N2C\\C1=C\\C)=Nc1c(OC)cc(OC)cc41 '
               'NAME: 19(E)-9,18-Didemethoxygardneramine REASON: WRONGLY '
               'CLASSIFIED Contains indole core, complex ring system (7 '
               'rings), and ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'COC(=O)C1=C2Nc3ccccc3[C@@]22CCN3CC[C@@H]4OC(=O)C[C@]4(C1)[C@@H]23 '
               'NAME: Apodine REASON: WRONGLY CLASSIFIED Contains indole core, '
               'complex ring system (6 rings), and ester group, ethyl/vinyl '
               'group, 3 sp2 carbons\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)[C@@H]4CCCN4C)C5=CC=CC=C5N2C)[C@@H](C)CO '
               'NAME: LSM-14052 REASON: WRONGLY CLASSIFIED Contains indole '
               'core, complex ring system (5 rings), and ethyl/vinyl group, 2 '
               'sp2 carbons\n'
               ' * SMILES: '
               'CC1=CC=C(C=C1)S(=O)(=O)N2CC3(C2)CN([C@@H](C4=C3C5=C(N4C)C=C(C=C5)OC)CO)CC6=CC=CC=C6 '
               'NAME: '
               "[(1S)-7-methoxy-9-methyl-1'-(4-methylphenyl)sulfonyl-2-(phenylmethyl)-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (6 rings), and ethyl/vinyl group\n'
               ' * SMILES: '
               'O=C1N2[C@@]3(C(=O)N[C@]14[C@@H](O[C@@]5(NC=6C7=C(OC(C)(C)C=C7)C=CC6C5=O)C([C@@H]4C3)(C)C)O)CCC2 '
               'NAME: Notoamide O REASON: WRONGLY CLASSIFIED Contains indole '
               'core, complex ring system (8 rings), and ethyl/vinyl group, 4 '
               'sp2 carbons\n'
               ' * SMILES: '
               'COCC(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3)C=C(C=C4)OC)CO)CC5=CC=CC=C5F '
               'NAME: '
               "1-[(1R)-2-[(2-fluorophenyl)methyl]-1-(hydroxymethyl)-7-methoxy-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]yl]-2-methoxyethanone "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (5 rings), and ethyl/vinyl group, 1 sp2 carbons\n'
               ' * SMILES: '
               'O=[N+]([O-])[C@@]12[C@@H](C(C)(C)[C@]3(C1)C(=O)N4C(CC(C=5C4=C3C=CC5OC)=O)(C)C)C[C@@H]6CCCN6C2 '
               'NAME: Clopiamine C REASON: WRONGLY CLASSIFIED Contains indole '
               'core, complex ring system (6 rings), and ethyl/vinyl group, 2 '
               'sp2 carbons\n'
               ' * SMILES: '
               'COC(=O)[C@@H]1C[C@]23CCCN2CC[C@]2(C(=O)OC)c4ccccc4N(C(=O)OC)[C@]12CC3 '
               'NAME: Grandilodine A, (rel)- REASON: WRONGLY CLASSIFIED '
               'Contains indole core, complex ring system (6 rings), and ester '
               'group, ethyl/vinyl group, 3 sp2 carbons\n'
               ' * SMILES: '
               'C=1(C[C@]2([C@@]3([NH+](C[C@H]4[C@@H]2O4)CC[C@]53C6=C(NC51)C=CC=C6)[H])[C@@H](C)OC(=O)C)C(=O)OC '
               'NAME: 19-O-acetylhoerhammericine(1+) REASON: WRONGLY '
               'CLASSIFIED Contains indole core, complex ring system (6 '
               'rings), and ester group, ethyl/vinyl group, 3 sp2 carbons\n'
               ' * SMILES: '
               'C1CCC(=CC1)C#CC2=CC3=C(C=C2)NC(=O)C34C(C5C(=O)OC(C(N5C4C6=CC=CC=C6OCCO)C7=CC=CC=C7)C8=CC=CC=C8)C(=O)NCC(C9=CC=CC=C9)O '
               'NAME: '
               "5-[2-(1-cyclohexenyl)ethynyl]-6'-[2-(2-hydroxyethoxy)phenyl]-N-(2-hydroxy-2-phenylethyl)-1',2-dioxo-3',4'-diphenyl-8'-spiro[1H-indole-3,7'-4,6,8,8a-tetrahydro-3H-pyrrolo[2,1-c][1,4]oxazine]carboxamide "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (9 rings), and ester group, ethyl/vinyl group, 4 sp2 '
               'carbons\n'
               ' * SMILES: '
               'O=C1N[C@@H](CC2=CNC=3C2=C(N[C@H]1C(C)C)C=C4[C@@](C(C)C)(CC[C@@](C34)(C=C)C)C)CO '
               'NAME: N-Demethylteleocidin B-4 REASON: WRONGLY CLASSIFIED '
               'Contains indole core, complex ring system (4 rings), and '
               'ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               '[C@@]123C(=NC4=C1C=CC(=C4)OC)[C@](C[C@]5([C@@]2([NH+](CC=C5)CC3)[H])CC)(C(=O)OC)O '
               'NAME: '
               '(3R)-1,2-didehydro-3-hydroxy-16-methoxy-2,3-dihydrotabersoninium '
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (5 rings), and ester group, ethyl/vinyl group, 3 sp2 '
               'carbons\n'
               ' * SMILES: '
               'CCCC(=O)N1CCC2(CC1)CN([C@@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)C(=O)C '
               'NAME: '
               "1-[(1S)-2-acetyl-1-(hydroxymethyl)-7-methoxy-9-methyl-1'-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]-1-butanone "
               'REASON: WRONGLY CLASSIFIED Contains indole core, complex ring '
               'system (4 rings), and ethyl/vinyl group, 2 sp2 carbons\n'
               ' * SMILES: '
               'C[C@H]1C[C@@H](C)\\C=C(CO)\\C[C@H](C)C(=O)N[C@@H](C)C(=O)N(C)[C@H](Cc2c(Br)[nH]c3ccccc23)C(=O)N[C@H](CC(=O)O1)c1ccc(O)cc1 '
               'NAME: jaspamide L REASON: WRONGLY CLASSIFIED Contains indole '
               'core, complex ring system (4 rings), and ester group, '
               'ethyl/vinyl group, 5 sp2 carbons\n'
               ' * SMILES: '
               'O=C/1N[C@H](C(=O)N\\C1=C/C=2C3=C(C(=CC=C3)CC4OC4(CO)C)NC2C(C=C)(C)C)C '
               'NAME: Variecolorin N REASON: WRONGLY CLASSIFIED Contains '
               'indole core, complex ring system (4 rings), and ethyl/vinyl '
               'group, 4 sp2 carbons\n'
               ' * SMILES: '
               'O=C1N2[C@@]3(C(=O)N[C@@]14[C@@H](C(C)(C)[C@@]5(C4)C(=O)NC6=C5C=CC7=C6C=CC(O7)(C)C)C3)CCC2 '
               'NAME: Taichunamide E REASON: WRONGLY CLASSIFIED Contains '
               'indole core, complex ring system (8 rings), and ethyl/vinyl '
               'group, 4 sp2 carbons\n'
               'False negatives: SMILES: '
               'COC=1C=CC=2C=3CCN4C[C@@]5(CC[C@@]([C@@]([C@]5(C[C@]4(C3NC2C1)[H])[H])(C(=O)OC)[H])(O)[H])[H] '
               'NAME: 11-methoxyyohimbine REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               'CCOC(=O)OC1=C(C=C(C=C1OC)C(=O)O[C@@H]2C[C@@H]3CN4CCC5=C([C@H]4C[C@@H]3[C@@H]([C@H]2OC)C(=O)OC)NC6=C5C=CC(=C6)OC)OC '
               'NAME: Syrosingopine REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               'COC(=O)[C@H]1[C@H]2C[C@H]3C4=C(CCN3C[C@@H]2CCC1O)C5=CC=CC=C5N4 '
               'NAME: '
               '(1S,15R,19S,20S)-18-hydroxy-1,3,11,12,14,15,16,17,18,19,20,21-dodecahydroyohimban-19-carboxylic '
               'acid methyl ester REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               'C/C=C(/C=O)\\[C@@]1(CC2C3=C(CCN2C)C4=CC=CC=C4N3[C@](C1)(O)[H])[H] '
               'NAME: Akagerine REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               '[H][C@@]12CCCC[C@@]1([H])C[C@]1([H])N(CCc3c1[nH]c1ccccc31)C2 '
               'NAME: yohimban REASON: MISSED Insufficient nitrogen-containing '
               'rings\n'
               ' * SMILES: '
               'C=C[C@]1(CN2CC[C@]1(C[C@]2(C3=C(CCO)C4=CC=CC=C4N3)[H])[H])[H] '
               'NAME: Cinchonamine HCl REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               'C/C=C\\1/CN2CC[C@@]34C5=CC=CC=C5N6[C@]4([C@@](C=C([C@]7(C[C@@]89C%10=CC=CC=C%10N%11C=C%12[C@]/%13(C[C@]%14([C@@]%15(CCN%14C\\C%13=C\\C)C%16=CC=CC=C%16N(C=C([C@]/%17(C[C@@]8(N7C\\C%17=C\\C)[H])[H])[C@@]9%11[H])[C@@]%12%15[H])[H])[H])[H])C6=O)([C@]1(C[C@@]32[H])[H])[H])[H] '
               'NAME: Strychnohexamine REASON: MISSED Complex structure lacks '
               'typical MIA ring patterns\n'
               ' * SMILES: '
               'CN1C=2C=CC(CC=3C=CC4=C(C3)[C@@]56CCN7CCC[C@@]8(CC([C@@]5(CC8)N4C)C(OC)=O)[C@@H]67)=CC2C9%10CCN%11CCC[C@@]%12(CC(C91CC%12)C(OC)=O)[C@@H]%10%11 '
               'NAME: Pleiokomenine A REASON: MISSED Complex structure lacks '
               'typical MIA ring patterns\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@H](O)C[C@@]1([H])C[C@]1([H])N(CCc3c1[nH]c1ccccc31)C2 '
               'NAME: 17alpha-yohimbol REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               ' * SMILES: '
               'CO[C@H]1[C@@H](C[C@@H]2CN3CCc4c([nH]c5ccccc45)[C@H]3C[C@@H]2[C@@H]1C(=O)OC)OC(=O)c1cc(OC)c(OC)c(OC)c1 '
               'NAME: deserpidine REASON: MISSED Insufficient '
               'nitrogen-containing rings\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Molecular weight 199.1 outside '
                                               'typical MIA range'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for alkaloid'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Molecular weight 1437.5 '
                                               'outside typical MIA range'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for alkaloid'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Too few rings (2) for MIA '
                                               'structure'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Too few rings (0) for MIA '
                                               'structure'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No indole or modified indole '
                                               'core found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for alkaloid'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No indole or modified indole '
                                               'core found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Molecular weight 220.2 outside '
                                               'typical MIA range'}],
    'sample_false_negatives': [   {   'smiles': 'C/C=C\\1/CN2CC[C@@]34C5=CC=CC=C5N6[C@]4([C@@](C=C([C@]7(C[C@@]89C%10=CC=CC=C%10N%11C=C%12[C@]/%13(C[C@]%14([C@@]%15(CCN%14C\\C%13=C\\C)C%16=CC=CC=C%16N(C=C([C@]/%17(C[C@@]8(N7C\\C%17=C\\C)[H])[H])[C@@]9%11[H])[C@@]%12%15[H])[H])[H])[H])C6=O)([C@]1(C[C@@]32[H])[H])[H])[H]',
                                      'name': 'Strychnohexamine',
                                      'reason': 'Complex structure lacks '
                                                'characteristic MIA patterns'},
                                  {   'smiles': 'CN1C=2C=CC(CC=3C=CC4=C(C3)[C@@]56CCN7CCC[C@@]8(CC([C@@]5(CC8)N4C)C(OC)=O)[C@@H]67)=CC2C9%10CCN%11CCC[C@@]%12(CC(C91CC%12)C(OC)=O)[C@@H]%10%11',
                                      'name': 'Pleiokomenine A',
                                      'reason': 'Complex structure lacks '
                                                'characteristic MIA patterns'},
                                  {   'smiles': 'CC[C@@]12CCCN3C(C=C(C4=CC=CC=C4NC(CC1)=O)[C@@]23O)=O',
                                      'name': 'Leuconolam',
                                      'reason': 'No indole or modified indole '
                                                'core found'},
                                  {   'smiles': 'C(C)[C@@]12[C@]3([C@]4(C=5C(N(C)[C@]4([C@](C(NNC(OCCSSC[C@H](NC([C@@H](NC([C@@H](NC([C@@H](NC([C@@H](NC(CC[C@H](NC(=O)C6=CC=C(NCC=7N=C8C(=NC(N)=NC8=O)NC7)C=C6)C(O)=O)=O)CC(O)=O)=O)CCCNC(=N)N)=O)CC(O)=O)=O)CC(O)=O)=O)C(O)=O)=O)=O)(O)[C@@H]1O)[H])=CC(OC)=C([C@]9(C(OC)=O)C%10=C(C=%11C(=CC=CC%11)N%10)CCN%12C[C@](C[C@](C%12)(CC)O)([H])C9)C5)CCN3CC=C2)[H]',
                                      'name': 'vintafolide',
                                      'reason': 'Molecular weight 1915.7 '
                                                'outside typical MIA range'},
                                  {   'smiles': 'C/C=C\\1/CN2[C@](C[C@@]34C5=CC=CC=C5N6C(C=CC([C@]1(C[C@@]32[H])[H])[C@@]46[H])=O)([C@]7([C@]8([C@@]9([C@]%10(C[C@]%11([C@@]%12(CCN%11CC%10=CCO8)C%13=CC=CC=C%13N([C@@]9%12[H])C7=O)[H])[H])[H])[H])[H])[H]',
                                      'name': 'Strychnogucine C',
                                      'reason': 'Complex structure lacks '
                                                'characteristic MIA patterns'},
                                  {   'smiles': '[H][C@@]12C[C@]3(CN4CC[C@@]56C(Nc7ccccc57)=C(C[C@@]5(CCO[C@@]35[H])[C@]46[H])C(=O)OC)CN3c4c(OC)cccc4[C@@]4(CCN5CC[C@]6([H])OCC[C@]6(C1)[C@@]45[H])[C@]23O',
                                      'name': 'vobtusine',
                                      'reason': 'Complex structure lacks '
                                                'characteristic MIA patterns'},
                                  {   'smiles': 'C1C[C@]23CC(C4(CC2)C5(CCN(C1)[C@@H]35)C6=CC(=CC=C6N4C)CC7=CC=C8C(=C7)[C@@]9%10CCN%11CCC[C@@]%12(CC([C@@]9(CC%12)N8)C(O)=O)[C@@H]%10%11)C(OC)=O',
                                      'name': 'Pleiokomenine B',
                                      'reason': 'Complex structure lacks '
                                                'characteristic MIA patterns'}],
    'attempt': 3,
    'reasoning': 'Based on the outcomes analysis, here are the key issues and '
                 'proposed improvements:\n'
                 '\n'
                 '1. False Negatives Analysis:\n'
                 '- Many legitimate MIAs were rejected due to "Insufficient '
                 'nitrogen-containing rings"\n'
                 '- The requirement for 2 nitrogen-containing rings is too '
                 'strict\n'
                 '- Some complex structures like Strychnohexamine were '
                 'rejected for lacking "typical MIA patterns"\n'
                 '- The pattern matching for nitrogen-containing rings needs '
                 'to be more flexible\n'
                 '\n'
                 '2. False Positives Analysis:\n'
                 '- Many compounds with indole cores were wrongly classified\n'
                 "- The current structural patterns aren't specific enough to "
                 'MIAs\n'
                 "- The complexity score isn't effectively distinguishing MIAs "
                 'from other complex indole-containing molecules\n'
                 '\n'
                 "Here's how we can improve the classifier:\n"
                 '\n'
                 '1. Revise Core Requirements:\n'
                 '- Relax the nitrogen ring requirement but make other '
                 'structural patterns more specific\n'
                 '- Add more specific MIA scaffold patterns based on the true '
                 'positives\n'
                 '- Include patterns for the secologanin-derived portion\n'
                 '\n'
                 '2. Improve Pattern Matching:\n'
                 '- Add more specific patterns for the terpene portion\n'
                 '- Include patterns for common MIA ring fusion types\n'
                 '- Add patterns for characteristic bridged ring systems\n'
                 '\n'
                 '3. Adjust Scoring:\n'
                 '- Increase weight for specific MIA patterns\n'
                 '- Add penalties for non-MIA structural features\n'
                 '- Make complexity score thresholds more selective\n'
                 '\n'
                 '4. Add New Features:\n'
                 '- Check for characteristic bridge atoms\n'
                 '- Look for specific ring fusion patterns common in MIAs\n'
                 '- Add patterns for common substituent positions\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 184,
    'num_false_positives': 100,
    'num_true_negatives': 4331,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.647887323943662,
    'recall': 0.9633507853403142,
    'f1': 0.7747368421052632,
    'accuracy': 0.9768498485504111,
    'negative_predictive_value': 0.9983863531581374}