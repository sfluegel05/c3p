"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_flavonoid, reason)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check carbon count (C15/C16 core + possible substituents)
    if c_count < 12:  # Relaxed from 15 to catch simpler flavonoids
        return False, "Too few carbons for flavonoid core structure"
        
    # Must contain oxygen
    if o_count < 1:  # Relaxed from 2
        return False, "Insufficient oxygen atoms for flavonoid"

    # Core structure patterns for different flavonoid classes
    patterns = [
        # Basic flavone/flavonol core (more flexible)
        "[#6]1~[#6]~[#6](~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2)~[#8]~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]12",
        # Isoflavone core
        "[#6]1~[#6]~[#6]2~[#8]~[#6]~[#6](=O)~[#6]2=C~[#6]1",
        # Chalcone pattern
        "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6](=O)~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1",
        # Flavanone core
        "[#6]1~[#6]~[#6](~[#6]2~[#6]~[#6]~[#6]~[#6]~[#6]2)~[#8]~[#6]2~[#6](=O)~[#6]~[#6]~[#6]12",
        # Anthocyanidin core
        "[#6]1~[#6]~[#6]2~[#8]~[#6](~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3)~[#6](~[#8])~[#6]~[#6]2~[#6]~[#6]1",
        # Aurone core
        "[#6]1~[#6]~[#6]2~[#8]~[#6](=C~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3)~[#6](=O)~[#6]2~[#6]~[#6]1",
        # More general diphenylpropane skeleton
        "[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1"
    ]
    
    # Check for presence of any core patterns
    found_core = False
    for pattern in patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            found_core = True
            break
            
    if not found_core:
        return False, "No flavonoid core structure found"

    # Look for typical oxygen-containing functional groups
    o_patterns = [
        # Ketone/carbonyl
        "[#6][#6](=[#8])[#6]",
        # Hydroxyl
        "[#8H1]",
        # Ether
        "[#6][#8][#6]",
        # Glycoside pattern
        "[#6][#8][#6]1[#8][#6]([#6][#8])[#6]([#8])[#6]([#8])[#6]1",
        # Pyranone ring
        "[#6]1[#6][#6](=[#8])[#8][#6]2[#6][#6][#6][#6][#6]12"
    ]
    
    oxygen_groups = 0
    for pattern in o_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None:
            matches = mol.GetSubstructMatches(pattern_mol)
            oxygen_groups += len(matches)
            
    if oxygen_groups == 0:  # Relaxed from 2
        return False, "No characteristic oxygen-containing functional groups"

    # Additional structural features common in flavonoids
    features = [
        # Phenol group
        "c1ccccc1[OH]",
        # Common substitution patterns
        "[#8]c1ccccc1",
        # Conjugated system
        "c1ccccc1",
        # Common glycosylation pattern
        "[#6][#8][#6]1[#8][#6][#6][#6][#6][#6]1[#8]"
    ]
    
    feature_count = 0
    for pattern in features:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            feature_count += 1

    if feature_count == 0:
        return False, "Missing typical flavonoid structural features"

    # Calculate number of rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Insufficient ring count for flavonoid structure"

    return True, "Contains flavonoid core structure with appropriate substitution patterns"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72544',
                          'name': 'flavonoids',
                          'definition': 'Any organic molecular entity whose '
                                        'stucture is based on derivatives of a '
                                        'phenyl-substituted 1-phenylpropane '
                                        'possessing a C15 or C16 skeleton, or '
                                        'such a structure which is condensed '
                                        'with a C6-C3 lignan precursors. The '
                                        "term is a 'superclass' comprising all "
                                        'members of the classes of flavonoid, '
                                        'isoflavonoid, neoflavonoid, '
                                        'chalcones, dihydrochalcones, aurones, '
                                        'pterocarpan, coumestans, rotenoid, '
                                        'flavonolignan, homoflavonoid and '
                                        'flavonoid oligomers. Originally '
                                        'restricted to natural products, the '
                                        'term is also applied to synthetic '
                                        'compounds related to them.',
                          'parents': ['CHEBI:25806', 'CHEBI:50860'],
                          'xrefs': ['Wikipedia:Flavonoids'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5327102803738317 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C(O)=C(C=C2)C(=O)/C=C/C3=CC(OC)=C(O)C=C3)CO '
               "NAME: 2',4',4-Trihydroxy-3',3-dimethoxychalcone 4'-O-glucoside "
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccccc1 NAME: '
               'pinocembrin REASON: CORRECT Contains flavonoid core structure '
               'with appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC=C(O)C=C4)COC(=O)C '
               'NAME: Kaempferol 3-(6-acetylgalactoside) REASON: CORRECT '
               'Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)CC(O)=O)O[C@@H](Oc2ccc3[C@@H]4Oc5cc6OCOc6cc5[C@@H]4COc3c2)[C@H](O)[C@H]1O '
               "NAME: (-)-Maackiain-3-O-glucosyl-6''-O-malonate REASON: "
               'CORRECT Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: '
               'O1C([C@H](O)[C@H](O)C(O)[C@@H]1OC=2C=C3OC(CC(=O)C3=C(O)C2)C4=CC=CC=C4)C '
               'NAME: Pinocembrin 7-O-rhamnoside REASON: CORRECT Contains '
               'flavonoid core structure with appropriate substitution '
               'patterns\n'
               ' * SMILES: '
               'O(C1C(O)[C@H](O)C(O[C@H]1OC=2C=C3C(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)COC(=O)CC(O)=O)=CC(O)=CC3=[O+]C2C5=CC=C(O)C=C5)COC(=O)/C=C/C6=CC=C(O)C=C6)[C@@H]7OC[C@@H](O)C(O)[C@H]7O '
               'NAME: Pelargonidin '
               "3-(6''-p-coumarylsambubioside)-5-(6'''-malonylglucoside) "
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1C(O)=O)C(C(=O)C2=C(OC)C3=C(OC=C3)C=C2)C(=O)C4=CC(O)=C(O)C=C4 '
               'NAME: '
               '6-[1-(3,4-dihydroxyphenyl)-3-(4-methoxy-1-benzofuran-5-yl)-1,3-dioxopropan-2-yl]-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=C(O)C=C(C=C2)C=3OC=4C(C(=O)C3OC)=C(O)C=5OCOC5C4)C(O)=O '
               "NAME: 3',4',5-Trihydroxy-3-methoxy-6,7-methylenedioxyflavone "
               "4'-glucuronide REASON: CORRECT Contains flavonoid core "
               'structure with appropriate substitution patterns\n'
               ' * SMILES: '
               'COc1ccc(cc1OC)[C@@H]1CC(=O)c2cc(CC=C(C)C)c3OC(C)(C)C=Cc3c2O1 '
               'NAME: '
               "(2S)-6-(gamma,gamma-dimethylallyl)-3',4'-dimethoxy-6'',6''-dimethylpyran[2'',3'':7,8]flavanone "
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(OC)C=C4)COC(=O)CC5CCNCC5 '
               'NAME: Isoflavone base + 1o, 1meo, o-hex+c7h12no REASON: '
               'CORRECT Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC=2C=C(C=CC2)C(=O)/C=C/C3=CC=CC=C3)C(O)=O '
               'NAME: '
               '3,4,5-trihydroxy-6-{3-[(2E)-3-phenylprop-2-enoyl]phenoxy}oxane-2-carboxylic '
               'acid REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1[C@H]([C@H](O)[C@@H](O)[C@@H](OC(=O)C=2C(OC)=CC=CC2)[C@@H]1OC=3C=C4OC(=C(O)C(=O)C4=C(O)C3)C5=CC=C(O)C=C5)C '
               'NAME: Gladiatoside C2 REASON: CORRECT Contains flavonoid core '
               'structure with appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1CO)C=2C(O)=C(C3OC(C(O)C(O)C3O)CO)C(O)=C4C2OC(=CC4=O)C=5C(OC6OC(C(O)C(O)C6O)C(O)=O)=CC(O)=C(O)C5 '
               'NAME: '
               '6-(2-{5,7-dihydroxy-4-oxo-6,8-bis[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]-4H-chromen-2-yl}-4,5-dihydroxyphenoxy)-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'C=1(C(C[C@@H](CC=C(C)C)C(C)=C)=C2O[C@H](C=3C(=CC=CC3O)O)CC(C2=C(O)C1)=O)OC '
               'NAME: exiguaflavanone B REASON: CORRECT Contains flavonoid '
               'core structure with appropriate substitution patterns\n'
               ' * SMILES: '
               'COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O[C@@H]1O[C@H](COC(=O)\\C=C\\c2ccc(O)cc2)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: malvidin 3-O-(6-O-(E)-4-coumaroyl-beta-D-glucoside) '
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=C(C=C2)C(=O)\\C=C\\C3=C(OC)C4=C(OC(C=C4)(C)C)C=C3)C(O)=O '
               'NAME: '
               '3,4,5-trihydroxy-6-{4-[3-(5-methoxy-2,2-dimethyl-2H-chromen-6-yl)prop-2-enoyl]phenoxy}oxane-2-carboxylic '
               'acid REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O[C@@H]1O[C@H](COC(=O)\\C=C/c2ccc(O)cc2)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: malvidin 3-O-(6-O-(Z)-4-coumaroyl-beta-D-glucoside) '
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(OC2OCC(O)C(O)C2O)C1OC=3C(=[O+]C=4C(C3)=C(O)C=C(O)C4)C5=CC(O)=C(O)C=C5)COC6OC(C(O)C(O)C6O)COC(=O)C7=CC=C(O)C=C7 '
               'NAME: Cyanidin '
               '3-O-[b-D-Xylopyranosyl-(1->2)-[(4-hydroxybenzoyl)-(->6)-b-D-glucopyranosyl-(1->6)]-b-D-galactopyranoside] '
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: O1CC(C2=C(OC)C=C(OC)C=C2)C(=O)C3=C1C=C(O)C=C3 NAME: '
               'Sativanone REASON: CORRECT Contains flavonoid core structure '
               'with appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C4=C(OC)C=5CCC(OC5C=C4O)(C)C)C(=O)C3C=C2)C(O)=O '
               'NAME: '
               '3,4,5-trihydroxy-6-{[3-(7-hydroxy-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-6-yl)-4-oxo-4H-chromen-7-yl]oxy}oxane-2-carboxylic '
               'acid REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(OC(=O)C2=CC=C(O)C=C2)C1OC3=CC=4OC=C(C(=O)C4C=C3)C5=CC=C(OC)C=C5)CO '
               'NAME: Formononetin 7-(2-p-hydroxybenzoylglucoside) REASON: '
               'CORRECT Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: '
               'O1C([C@H](OC(=O)C)C(OC(=O)C)C(OC(=O)C)[C@@H]1OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)CO '
               "NAME: Quercetin 3-(2'',3'',4''-triacetylgalactoside) REASON: "
               'CORRECT Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: COc1ccc(cc1)[C@@H]1CC(=O)c2c(O)c(C)c(O)c(C)c2O1 '
               'NAME: matteucinol REASON: CORRECT Contains flavonoid core '
               'structure with appropriate substitution patterns\n'
               ' * SMILES: '
               'O1[C@H](C(O)=O)[C@H]([C@H](O)[C@H]([C@@H]1OC=2C=C(C3=C(OC(CC3=O)C=4C=CC(=CC4)O)C2)O)O)O '
               'NAME: naringenin 7-O-beta-D-glucuronide REASON: CORRECT '
               'Contains flavonoid core structure with appropriate '
               'substitution patterns\n'
               ' * SMILES: '
               'O1C(CC(=O)C=2C1=C(CC=C(C)C)C(O)=CC2O)C3=CC(=C(O)C(O)=C3)CC=C(C)C '
               "NAME: (S)-3',4',5,7-Tetrahydroxy-5',8-diprenylflavanone "
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substitution patterns\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O(C1C(O)[C@H](O)C(O[C@H]1OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(OC)=C(O)C=C4)CO)[C@@H]5OC[C@@H](O)C(O)[C@H]5O '
               'NAME: Isorhamnetin 3-sambubioside REASON: MISSED Insufficient '
               'oxygen-containing functional groups\n'
               ' * SMILES: C1=C2C(=C(C(=C1)O)O)OC(=CC2=O)C=3C=CC=CC3 NAME: '
               '7,8-dihydroxyflavone REASON: MISSED No flavonoid core '
               'structure found\n'
               ' * SMILES: '
               'CC1=CC(N2C3=C1C=CC=C3C4(C2=O)C5=C(C=C(C=C5)O)OC(=C4C#N)N)(C)C '
               'NAME: LSM-26505 REASON: MISSED No flavonoid core structure '
               'found\n'
               ' * SMILES: '
               'COc1ccc(cc1O)-c1coc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)ccc2c1=O '
               'NAME: calycosin-7-O-beta-D-glucoside REASON: MISSED '
               'Insufficient oxygen-containing functional groups\n'
               ' * SMILES: O1C2=C(C(=O)C(C3=CC=C(OC)C=C3)=C1)C=CC(OC)=C2O '
               'NAME: Retusin 7-Methyl Ether REASON: MISSED No flavonoid core '
               'structure found\n'
               ' * SMILES: '
               'O=C1C=2C=3C=4C(OC=5C=CC=C(C5C4[C@@H](O)[C@@H](C3C(=C1)O)C)O)=CC2O '
               'NAME: (-)-dalmanol A REASON: MISSED No flavonoid core '
               'structure found\n'
               ' * SMILES: O(C1=CC=C(C(=O)\\C=C\\C=2C(O)=CC=CC2)C=C1)C NAME: '
               '(2E)-3-(2-hydroxyphenyl)-1-(4-methoxyphenyl)prop-2-en-1-one '
               'REASON: MISSED No flavonoid core structure found\n'
               ' * SMILES: '
               'O(C1C(O)C(O)C(OC1C)OC=2C(=O)C=3C(OC2C4=CC(O)=C(O)C=C4)=CC(OC5OC(C(O)C(O)C5O)COC6OC(C(O)C(O)C6O)C)=CC3O)C7OC(C(O)C(O)C7O)CO '
               'NAME: Quercetin 3-(glucosyl-(1->4)-rhamnoside) 7-rutinoside '
               'REASON: MISSED Insufficient oxygen-containing functional '
               'groups\n'
               ' * SMILES: '
               'COc1cc(O)cc2oc(-c3cc(O)c(O)c(O)c3)c(O[C@@H]3O[C@@H](C)[C@H](O)[C@@H](O)[C@H]3O)c(=O)c12 '
               'NAME: myricitrin-5-methyl ether REASON: MISSED Insufficient '
               'oxygen-containing functional groups\n'
               ' * SMILES: '
               'Oc1cc(O)c2c(c1)oc(-c1ccc(O)c(O)c1)c(OS(O)(=O)=O)c2=O NAME: '
               'quercetin 3-sulfate REASON: MISSED No flavonoid core structure '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@H](COC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCCCCCCCC',
                                     'name': 'TG(14:0/14:0/22:2(13Z,16Z))',
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'O([C@@H]1O[C@@H]([C@H](O)[C@H](O[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H]1O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)CO)[C@H]4[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]4CO)OC[C@@H](O)[C@H](O)[C@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@@H](NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5R,6R)-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-[(2R,3S,4R,5S)-5-acetamido-2,3,6-trihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyhexoxy]-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'P12(O[C@@]3(NC(=O)C)[C@@](O1)(O)[C@H](O)[C@H](O[C@]3(O2)O)CO)=O',
                                     'name': 'N-Acetylglucosamine phosphate',
                                     'reason': 'Too few carbons for flavonoid '
                                               'core structure'},
                                 {   'smiles': 'OC1=C(N=CC=C1)C',
                                     'name': '3-hydroxyl-2-methylpyridine',
                                     'reason': 'Too few carbons for flavonoid '
                                               'core structure'},
                                 {   'smiles': '[H][C@@]1(C)O[C@@]1([H])C1=C[C@H](O)[C@@H](C)OC1=O',
                                     'name': 'aspyrone',
                                     'reason': 'Too few carbons for flavonoid '
                                               'core structure'},
                                 {   'smiles': 'O=C1C(O)=C(NC2=C1C=C(C)C(=C2)C)C(=O)O',
                                     'name': '3,4-dihydroxy-6,7-dimethyl-quinoline-2-carboxylic '
                                             'acid',
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'O1[C@@]2(OC(C)(C)[C@H]([C@H]2O)C)[C@H]([C@@H]3[C@@]4([C@@]1(C5=C(C(=C(O)C=C5)C)C=C4)CC3)C)C',
                                     'name': 'Blazeispirol D',
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'CCS(=O)(=O)N1CC2(C1)CN([C@@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)C',
                                     'name': "[(1S)-1'-ethylsulfonyl-7-methoxy-2,9-dimethyl-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol",
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'CN(C)CC(=O)N1CCC2(CC1)CN([C@@H](C3=C2C4=C(N3)C=C(C=C4)OC)CO)C(=O)C5=CN=CC=C5',
                                     'name': "2-(dimethylamino)-1-[(1S)-1-(hydroxymethyl)-7-methoxy-2-[oxo(3-pyridinyl)methyl]-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]ethanone",
                                     'reason': 'No flavonoid core structure '
                                               'found'},
                                 {   'smiles': 'O=C1C2(C(C(CC1(C(O)=C(C2=O)CC=C(C)C)CC=C(C)C)CC=C(C)C)(CCC=C(C)C)C)C(=O)C(CC)C',
                                     'name': 'Adhyperforin',
                                     'reason': 'No flavonoid core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC1=CC(N2C3=C1C=CC=C3C4(C2=O)C5=C(C=C(C=C5)O)OC(=C4C#N)N)(C)C',
                                      'name': 'LSM-26505',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1C2=C(C(C3=CC(OC)=C(OC)C=C3)=CC1=O)C(OC)=CC(OC)=C2O',
                                      'name': "8-Hydroxy-5,7,3',4'-tetramethoxy-4-phenylcoumarin",
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'C12=C(C3=C(C(=CC(O3)=O)C4=CC=CC=C4)C5=C1C=CC(O5)(C)C)C(C(C)C(O2)C)O',
                                      'name': 'inophyllolide chromanol',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'CC[C@H](C)C(=O)c1c(O)c(CC=C(C)C)c2oc(=O)cc(-c3ccccc3)c2c1O',
                                      'name': 'Mammea A/AB, (-)-',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C=CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)C(C)C)(CO)C',
                                      'name': '5-hydroxy-8-(hydroxymethyl)-8-methyl-6-(2-methylpropanoyl)-4-phenyl-2H,8H-pyrano[2,3-f]chromen-2-one',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'COc1cc(ccc1O)C1CCOc2cc(OC)c(O)cc12',
                                      'name': 'pterolinus J',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1C(C=CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)CCC)(C)C',
                                      'name': '6-Butyryl-5-hydroxy-4-phenylseselin',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'BrC1=C(O)C=CC(=C1)[C@H]2C3=C(OC[C@H]2OC4=CC(O)=CC=C4)C(=CC=C3O)[C@@H](O)CN',
                                      'name': 'Ammonificin A',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'C=1(C=C(C=CC1O)CC(CCC(OCOS(=O)(=O)O)=O)O)O',
                                      'name': '4-Hydroxy-5-(dihydroxyphenyl)-valeric '
                                              'acid-O-methyl-O-sulphate',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1C(CC(C2=C1C=C(O)C3=C2N(C4=C(C3=O)C=CC=C4O)C)C5=C6OC(C=CC6=C7N(C8=C(C(=O)C7=C5O)C=CC=C8O)C)(C)C)(C)C',
                                      'name': 'Bis(5-hydroxynoracronycine)',
                                      'reason': 'No flavonoid core structure '
                                                'found'}],
    'attempt': 2,
    'reasoning': 'Based on the outcomes analysis, here are the key issues with '
                 'the previous code:\n'
                 '\n'
                 '1. False Negatives Analysis:\n'
                 '- The code missed several legitimate flavonoids like '
                 'Isorhamnetin 3-sambubioside and '
                 'calycosin-7-O-beta-D-glucoside due to overly strict oxygen '
                 'functional group requirements\n'
                 '- Simple flavonoids like 7,8-dihydroxyflavone were missed '
                 'due to insufficient core structure patterns\n'
                 '- Complex glycosylated flavonoids were missed due to the '
                 'molecular weight and pattern matching being too restrictive\n'
                 '\n'
                 '2. Core Structure Detection Issues:\n'
                 '- The SMARTS patterns for core structures were too specific '
                 'and rigid\n'
                 '- Need more flexible patterns to catch variations in the '
                 'basic flavonoid skeleton\n'
                 '- Should include more comprehensive patterns for different '
                 'flavonoid subclasses\n'
                 '\n'
                 '3. Improvements Needed:\n'
                 '- Add more flexible core structure patterns\n'
                 '- Relax the oxygen functional group requirements\n'
                 '- Better handle glycosylated flavonoids\n'
                 '- Add patterns for specific flavonoid subclasses like:\n'
                 '  * Isoflavones\n'
                 '  * Chalcones\n'
                 '  * Flavanones\n'
                 '  * Anthocyanins\n'
                 '  * Aurones\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 3680,
    'num_false_positives': 100,
    'num_true_negatives': 2713,
    'num_false_negatives': 66,
    'num_negatives': None,
    'precision': 0.9735449735449735,
    'recall': 0.9823812066203951,
    'f1': 0.9779431304809991,
    'accuracy': 0.9746912639121817,
    'negative_predictive_value': 0.9762504498020871}