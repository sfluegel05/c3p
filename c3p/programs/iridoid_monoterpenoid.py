"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:24532 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoids typically have a cyclopentane ring fused to a six-membered oxygen heterocycle,
    or their seco-derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Molecular weight check (typical range for iridoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 800:
        return False, "Molecular weight outside typical range for iridoids"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 30:
        return False, "Carbon count outside typical range for iridoids"
    
    if o_count < 2:
        return False, "Insufficient oxygen atoms for iridoid structure"

    # Define various iridoid skeleton patterns
    patterns = [
        # Classic iridoid skeleton
        "[C]1[C][C][C]2[C]1[O][C][C][C][C]2",
        # Seco-iridoid skeleton
        "[C]-[C]-[C]-[C]1[O][C][C][C][C]1",
        # Modified iridoid skeleton with oxygen bridge
        "[C]1[C][C]2[O][C][C][C][C]2[O]1",
        # Rearranged iridoid skeleton
        "[C]1[C]2[C][C][C]1[O][C][C][C]2",
        # Cyclopentane fused to pyran (more general)
        "[C]1[C][C]2[C][C]1[O][C][C][C]2"
    ]

    is_iridoid_skeleton = False
    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            is_iridoid_skeleton = True
            break

    if not is_iridoid_skeleton:
        # Check for secoiridoid patterns
        seco_patterns = [
            # Common secoiridoid patterns
            "[C]-[C]-[C]-[C]-[C]-[O]-[C]",
            "[C]-[C](=O)-[C]-[C]-[O]-[C]",
            "[C]-[C](=O)-[C]-[C]-[C]-[O]"
        ]
        
        for pattern in seco_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat and mol.HasSubstructMatch(pat):
                is_iridoid_skeleton = True
                break

    if not is_iridoid_skeleton:
        return False, "Missing characteristic iridoid/secoiridoid skeleton"

    # Check for characteristic substituent patterns
    substituent_patterns = [
        (Chem.MolFromSmarts("[CH2][OH]"), "hydroxymethyl"),
        (Chem.MolFromSmarts("[C](=O)[OH]"), "carboxyl"),
        (Chem.MolFromSmarts("[OH]"), "hydroxyl"),
        (Chem.MolFromSmarts("[C](=O)[O][C]"), "ester"),
        (Chem.MolFromSmarts("[CH3]"), "methyl")
    ]
    
    found_groups = []
    for pattern, name in substituent_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            found_groups.append(name)
    
    if len(found_groups) < 2:
        return False, "Insufficient characteristic substituents"

    # Ring count check
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2 or ring_count > 6:
        return False, "Ring count outside typical range for iridoids"

    # Exclude compounds that are likely to be simple glycosides
    sugar_pattern = Chem.MolFromSmarts("[OH]C1[CH]([OH])[CH]([OH])[CH]([OH])[CH]([OH])[CH]1[OH]")
    if sugar_pattern and mol.HasSubstructMatch(sugar_pattern):
        # Additional check to ensure it's not just a simple glycoside
        if ring_count == 2 and o_count > 5:
            return False, "Likely a simple glycoside"

    # Success case
    structure_type = "secoiridoid" if not any("cyclopentane" in p for p in patterns) else "iridoid"
    substituents = ", ".join(found_groups)
    return True, f"Identified as {structure_type} with {substituents} groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50563',
                          'name': 'iridoid monoterpenoid',
                          'definition': 'One of a class of monoterpenoids '
                                        'biosynthesized from isoprene and '
                                        'often intermediates in the '
                                        'biosynthesis of alkaloids. Iridoids '
                                        'usually consist of a cyclopentane '
                                        'ring fused to a six-membered oxygen '
                                        'heterocycle; cleavage of a bond in '
                                        'the cyclopentane ring gives rise to '
                                        'the subclass known as secoiridoids.',
                          'parents': ['CHEBI:25409'],
                          'xrefs': ['Wikipedia:Iridoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.3076923076923077 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O[C@@]1([C@@]2([C@](CC1)(C(=CO[C@H]2O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@H]3O)CO[C@@H]4OC[C@@](O)([C@@H]4O)CO)C(O)=O)[H])[H])C '
               'NAME: 6-O-beta-D-Apiofuranosyl-mussaenosidic acid REASON: '
               'CORRECT Identified as seco-iridoid with hydroxyl, carboxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O1[C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@]3([C@]([C@H](O)[C@H](OC(=O)/C=C/C4=CC=CC=C4)[C@@H]3C)(C(=C1)C(OC)=O)[H])[H] '
               'NAME: methyl '
               '(1S,4aS,5S,6R,7R,7aR)-5-hydroxy-7-methyl-6-[(E)-3-phenylprop-2-enoyl]oxy-1-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,4a,5,6,7,7a-hexahydrocyclopenta[c]pyran-4-carboxylate '
               'REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: O1C(OC)C2C([C@H](OC)C[C@@]2(O)C)CC1OC NAME: '
               'Mioporosidegenin REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, methyl groups\n'
               ' * SMILES: '
               '[H][C@]12C[C@H](O)[C@H](C)[C@@]1([H])[C@@H](OC=C2C(=O)OC)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: loganin REASON: CORRECT Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O([C@@H]1OC=C([C@@]2([C@]1(C(=C[C@@H]2O)COC(=O)C)[H])[H])C(OC)=O)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: Daphylloside REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: O(C1OC=C(C2C1C(C(O)C2)C)C(O)=O)C3OC(C(O)C(O)C3O)CO '
               'NAME: '
               '6-hydroxy-7-methyl-1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,4a,5,6,7,7a-hexahydrocyclopenta[c]pyran-4-carboxylic '
               'acid REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'carboxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'OC1(C2C(CC1O)C(=COC2OC3OC(C(O)C(O)C3O)COC(=O)/C=C/C4=CC=C(O)C=C4)C(O)=O)C '
               'NAME: Lippioside I REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, carboxyl, methyl groups\n'
               ' * SMILES: '
               '[H][C@]12C=C[C@](O)(CO)[C@@]1([H])[C@@H](OC=C2C(O)=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: monotropein REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, carboxyl, hydroxymethyl groups\n'
               ' * SMILES: '
               '[H][C@@]12CCOC(=O)C1=CO[C@@H](O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1OC(=O)c1c(O)cc(O)cc1-c1cccc(O)c1)[C@@H]2C=C '
               'NAME: amarogentin REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, hydroxymethyl groups\n'
               ' * SMILES: OC12C(C(O)(CC1O)C)C(OC3OC(C(O)C(O)C3O)CO)OC=C2 '
               'NAME: '
               '(1R,4Ar,5S,7R)-7-methyl-1-[(2R,3S,4R,5R,6S)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,5,6,7a-tetrahydrocyclopenta[c]pyran-4a,5,7-triol '
               'REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: O(C1OC=C(C2C1C(=CC2)CO)C(OC)=O)C3OC(C(O)C(O)C3O)CO '
               'NAME: Methyl '
               '7-(hydroxymethyl)-1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,4a,5,7a-tetrahydrocyclopenta[c]pyran-4-carboxylate '
               'REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: O(C1OC=C(C2C1C(O)(CC2)C)C(OC)=O)C3OC(C(O)C(O)C3O)CO '
               'NAME: '
               '5-Chloro-2-(4-(2-(2-hydroxyethoxy)ethyl)piperazin-1-yl)-N-(6-(3-methoxybenzamido)pyridin-3-yl)isonicotinamide '
               'REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O([C@@H]1OC=C([C@@]2([C@]1(C(=C[C@@H]2O)COC(=O)C)[H])[H])C(O)=O)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: Asperulosidic acid REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, carboxyl, methyl, hydroxymethyl '
               'groups\n'
               ' * SMILES: '
               'O1[C@@H](O[C@@H]2OC=C([C@H]([C@H]2C=C)CC(O)=O)C(OC)=O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1CO '
               'NAME: secoxyloganin REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, carboxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'CO[C@H]1[C@H](O[C@@H]2OC=C([C@H]3[C@@H]4O[C@@H]4[C@]4(OC(=O)C(=C4)C(=O)c4ccc(OC)c(O)c4)[C@@H]23)C(=O)OC)O[C@H](CO)[C@@H](O)[C@@H]1O '
               'NAME: citrifolinoside REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O([C@@H]1OC=C[C@@]2([C@]1([C@@](O)(C[C@H]2O)C)[H])[H])[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: Ajugol REASON: CORRECT Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O1[C@@H]2C=C(CC[C@]23[C@@H](CC[C@H]3C)[C@](C1)(O)C)C NAME: '
               'Cordycepol A REASON: CORRECT Identified as seco-iridoid with '
               'hydroxyl, methyl groups\n'
               ' * SMILES: '
               'OC12C(C(O)(CC1O)C)C(OC=C2)OC3OC(C(O)C(O)C3O)COC4OC(C(O)C(O)C4O)CO '
               'NAME: Stachyoside A REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               '[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)O[C@H]2[C@@H]([C@@H](C(=CO2)C(OC)=O)CC=O)C=C '
               'NAME: (-)-secologanin REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: OC1(C2C(CC1O)C(=COC2OC3OC(C(O)C(O)C3O)CO)C(O)=O)C '
               'NAME: Caryoptosidic acid REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, carboxyl, methyl, hydroxymethyl '
               'groups\n'
               ' * SMILES: '
               'O([C@@]1([C@@]2([C@](O)([C@H](O)C1)C=CO[C@H]2O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)[H])C)C(=O)C '
               'NAME: 8-O-Acetylharpagide REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O([C@]12[C@]([C@H](O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)OC=C1)(C(=C[C@H]2O)CO)[H])[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO '
               'NAME: Melittoside REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O[C@@]1([C@@]2([C@]([C@H](O)C1)(C(=CO[C@H]2O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO)C(OC)=O)[H])[H])C '
               'NAME: Shanzhiside methyl ester REASON: CORRECT Identified as '
               'seco-iridoid with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=CC3)O)=C/C '
               'NAME: ligstroside REASON: CORRECT Identified as seco-iridoid '
               'with hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O(C1OC=C(C2C1C(CC2)C)C(O)=O)[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO '
               'NAME: '
               '7-Methyl-1-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,4a,5,6,7,7a-hexahydrocyclopenta[c]pyran-4-carboxylic '
               'acid REASON: CORRECT Identified as seco-iridoid with hydroxyl, '
               'carboxyl, methyl, hydroxymethyl groups\n'
               'False positives: SMILES: '
               'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Identified as seco-iridoid '
               'with hydroxyl, carboxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N '
               'NAME: LSM-4497 REASON: WRONGLY CLASSIFIED Identified as '
               'seco-iridoid with hydroxyl, methyl groups\n'
               ' * SMILES: '
               'C1CCN(CC1)CCNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@H](O2)CO)OC4=C3C=C(C=C4)NC(=O)NC5CCCC5 '
               'NAME: '
               '2-[(1R,3S,4aS,9aR)-6-[[(cyclopentylamino)-oxomethyl]amino]-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[2-(1-piperidinyl)ethyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O '
               'NAME: beta-D-Galp-(1->2)-beta-D-Galp REASON: WRONGLY '
               'CLASSIFIED Identified as seco-iridoid with hydroxyl, '
               'hydroxymethyl groups\n'
               ' * SMILES: '
               'CNC(=O)C[C@H]1CC[C@@H]2[C@@H](O1)COC[C@H](CN2C(=O)NC3=CC(=CC=C3)OC)O '
               'NAME: '
               '(3S,6aR,8R,10aR)-3-hydroxy-N-(3-methoxyphenyl)-8-[2-(methylamino)-2-oxoethyl]-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocine-1-carboxamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl groups\n'
               ' * SMILES: '
               'CN(C)C1=CC2=C(C=C1)O[C@H]3[C@@H]2C[C@H](O[C@H]3CO)CC(=O)NCC4=CC=CC=C4OC '
               'NAME: '
               '2-[(1S,3S,4aR,9aS)-6-(dimethylamino)-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[(2-methoxyphenyl)methyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               '[C@@H]1(O)[C@H](O)[C@H]([C@H](O[C@@H]2[C@H]([C@@H](O[C@@H]([C@H]2O[C@@H]3O[C@H]([C@H]([C@H]([C@@H]3O)O)O)C)CO)OCCCCCCCCC(OC)=O)NC(C)=O)O[C@@H]1CO)O '
               'NAME: methyl '
               '8-{[beta-D-galactosyl-(1->3)-[alpha-L-fucosyl-(1->4)]-N-acetyl-beta-D-glucosaminyl]oxy}nonanoate '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O=C1C=C2[C@]3(OOC([C@@H]([C@@H]2O)O3)(C)C)[C@@H]4[C@]1(O4)C '
               'NAME: Acremine P REASON: WRONGLY CLASSIFIED Identified as '
               'seco-iridoid with hydroxyl, methyl groups\n'
               ' * SMILES: '
               'C1C[C@@H]2[C@H](COC[C@@H](CN2CC3=CC=C(C=C3)OC4=CC=CC=C4)O)O[C@H]1CC(=O)N5CCC6=CC=CC=C6C5 '
               'NAME: '
               '2-[(3R,6aR,8R,10aR)-3-hydroxy-1-[(4-phenoxyphenyl)methyl]-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocin-8-yl]-1-(3,4-dihydro-1H-isoquinolin-2-yl)ethanone '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl groups\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](O)[C@@H](O[C@@H]2CO)O)[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO '
               'NAME: '
               '(2S,3R,4S,5S,6R)-2-[[(2R,3S,4S,5R,6R)-6-[(2R,3S,4R,5R,6R)-4,5-Dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-6-(hydroxymethyl)oxane-3,4,5-triol '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'CN(C)C1=CC2=C(C=C1)O[C@H]3[C@@H]2C[C@@H](O[C@@H]3CO)CC(=O)NCC4=CC=CC=N4 '
               'NAME: '
               '2-[(1R,3R,4aR,9aS)-6-(dimethylamino)-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-(2-pyridinylmethyl)acetamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'CC(C)=CCNc1ncnc2n(cnc12)[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: 9-(alpha-D-glucosyl)-N(6)-isopentenyladenine REASON: '
               'WRONGLY CLASSIFIED Identified as seco-iridoid with hydroxyl, '
               'methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC=2C(=C(O)C=CC2)C(OCC3=CC=CC=C3)=O)CO NAME: '
               'Benzyl 2,6-dihydroxybenzoate 2-glucoside REASON: WRONGLY '
               'CLASSIFIED Identified as seco-iridoid with hydroxyl, '
               'hydroxymethyl groups\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC(=O)C2=CC=CC=C2)CO '
               'NAME: benzoyl-beta-D-glucoside REASON: WRONGLY CLASSIFIED '
               'Identified as seco-iridoid with hydroxyl, hydroxymethyl '
               'groups\n'
               ' * SMILES: '
               'O1C23C(O)(C(C4(O)C2(O)C(=C(C4)C(C)C)C)(CC1=O)C)CCC(C3O)C NAME: '
               'Anhydrocinnzeylanol REASON: WRONGLY CLASSIFIED Identified as '
               'seco-iridoid with hydroxyl, methyl groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C=C(OC)C(=C2[C@H](OC)[C@@H]3[C@H]1CO[C@](C3)(OC)C)O '
               'NAME: 3,5-dimethoxydihydrofusarubin B REASON: WRONGLY '
               'CLASSIFIED Identified as seco-iridoid with hydroxyl, methyl '
               'groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C3=C(C[C@](O)(C)C[C@@H]3O[C@@H]4O[C@H]([C@H](O)[C@H]([C@H]4O)O)C)C(=C2C(=O)C=5C1=C(O)C(=CC5)C)O '
               "NAME: 4-O,3'-O-didemethyl mutactimycin REASON: WRONGLY "
               'CLASSIFIED Identified as seco-iridoid with hydroxyl, methyl '
               'groups\n'
               ' * SMILES: '
               'COC1=CC=CC=C1CNC(=O)C[C@H]2C[C@@H]3[C@H]([C@H](O2)CO)OC4=C3C=C(C=C4)NS(=O)(=O)C5=CC=C(C=C5)F '
               'NAME: '
               '2-[(1R,3R,4aS,9aR)-6-[(4-fluorophenyl)sulfonylamino]-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[(2-methoxyphenyl)methyl]acetamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'C1CCC(CC1)NC(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@@H](O[C@H]4CO)CC(=O)NCC5=CN=CC=C5 '
               'NAME: '
               '2-[(1S,3R,4aS,9aR)-6-[[(cyclohexylamino)-oxomethyl]amino]-1-(hydroxymethyl)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-(3-pyridinylmethyl)acetamide '
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O(C1C(O)C(O)C(OC1OC)COC(=O)CCC(OC)=O)[C@@H]2[C@H](OC3=C2C=CC(O)=C3)CC4=CC=C(O)C=C4 '
               'NAME: Ambofuracin REASON: WRONGLY CLASSIFIED Identified as '
               'seco-iridoid with hydroxyl, methyl groups\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1OC(C(O)C)C2=CC=C(OC)C=C2)CO NAME: '
               '(-)-erythro-Anethole glycol 1-glucoside REASON: WRONGLY '
               'CLASSIFIED Identified as seco-iridoid with hydroxyl, methyl, '
               'hydroxymethyl groups\n'
               ' * SMILES: '
               'COC1=CC2=C(C=C1)C3=C(N2)[C@@H](N(CC34CCN(CC4)CC5CCOCC5)C(=O)NC6=CC=C(C=C6)F)CO '
               'NAME: '
               "(1R)-N-(4-fluorophenyl)-1-(hydroxymethyl)-7-methoxy-1'-(4-oxanylmethyl)-2-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]carboxamide "
               'REASON: WRONGLY CLASSIFIED Identified as seco-iridoid with '
               'hydroxyl, methyl, hydroxymethyl groups\n'
               ' * SMILES: '
               'O1[C@]([C@]2([C@]3(OC(=O)C)[C@](CC2)(CC=C(CC3)C(OC)=O)C1=O)[H])(C)/C=C/C=C(\\C)/C(O)=O '
               'NAME: Pseudolaric acid B REASON: WRONGLY CLASSIFIED Identified '
               'as seco-iridoid with hydroxyl, carboxyl, methyl groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C(=CC=C2C(=O)C3=C1C=C(O)C=4C3=C(O)C=C(C)C4)[C@@H]5O[C@H]([C@@H](N(C)C)CC5)C '
               'NAME: Monacyclinone B REASON: WRONGLY CLASSIFIED Identified as '
               'seco-iridoid with hydroxyl, methyl groups\n'
               ' * SMILES: '
               'O1[C@H]([C@H](NC(=O)CO)[C@@H](O)C[C@]1(OC[C@H]2OC(O)[C@H](NC(=O)C)[C@@H](OC3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H]2O)C(O)=O)[C@H](O)[C@H](O)CO '
               'NAME: '
               '(2R,4S,5R,6R)-2-[[(2R,3R,4R,5R)-5-Acetamido-3,6-dihydroxy-4-[(3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]methoxy]-4-hydroxy-5-[(2-hydroxyacetyl)amino]-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Identified as seco-iridoid '
               'with hydroxyl, carboxyl, methyl, hydroxymethyl groups\n'
               'False negatives: SMILES: OC12C(C(CC1(CO)C)(C)C)CCC(C2C(O)=O)C '
               'NAME: '
               '3a-hydroxy-3-(hydroxymethyl)-1,1,3,5-tetramethyl-2,4,5,6,7,7a-hexahydroindene-4-carboxylic '
               'acid REASON: MISSED Missing required pyran or modified oxygen '
               'heterocycle\n'
               ' * SMILES: '
               'ClCC1(O)C2C(=CC1OC(=O)CC(C)C)C(=COC2OC(=O)CC(C)C)COC(=O)C '
               'NAME: Valechlorin REASON: MISSED Missing required pyran or '
               'modified oxygen heterocycle\n'
               ' * SMILES: COC(=O)C1=CO[C@@H](O)[C@H]2[C@@H]1CC=C2CO NAME: '
               'Genipin REASON: MISSED Missing required pyran or modified '
               'oxygen heterocycle\n'
               ' * SMILES: '
               'C12=C3[C@H]([C@@H]4C(CC[C@@]4(C3)C)(C)C)CC[C@]1([C@H](C)CC2)C '
               'NAME: Cattleyene REASON: MISSED No oxygen atoms found\n'
               ' * SMILES: '
               'C[C@@H]1C[C@@]2(CC(C)=CCC[C@@H]2[C@@]1(CO)CC\\C(CO)=C\\CO)C=O '
               'NAME: Portulal REASON: MISSED Missing required pyran or '
               'modified oxygen heterocycle\n'
               ' * SMILES: '
               'O=C1NCCC[C@@H]2NC(=O)C(C2=O)=C(O)C=C[C@H]3[C@H](CC=C1)[C@H]4[C@H](O)[C@@H](O)[C@H]5[C@H]([C@H](O)C)[C@@H](C[C@@H]5[C@H]4C3)C '
               'NAME: Hydroxyikarugamycin C REASON: MISSED Missing required '
               'pyran or modified oxygen heterocycle\n'
               ' * SMILES: C1(=C[C@H]2[C@@]3([C@H](CC[C@@H]3C2(C)C)C)CC1)C '
               'NAME: Italicene REASON: MISSED No oxygen atoms found\n'
               ' * SMILES: '
               'O1C2(C3C(=CC2OC(=O)CC(C)C)C(=COC3OC(=O)CC(C)C)COC(=O)CC(O)(C)C)C1 '
               'NAME: 11-Deacetylvaltrate 11-(3-hydroxy-3-methylbutanoate) '
               'REASON: MISSED Missing required pyran or modified oxygen '
               'heterocycle\n'
               ' * SMILES: '
               '[H][C@]12C[C@H](O)[C@H](C)[C@@]1([H])[C@H](O)OC=C2C(=O)OC '
               'NAME: loganetin REASON: MISSED Missing required pyran or '
               'modified oxygen heterocycle\n'
               ' * SMILES: '
               'Cl/C=C\\1/C2=C3C(=C[C@]4(C[C@@H]([C@@]([C@H]4[C@@H]3CC1)(C[C@H](O)[C@H](O)[C@@](O)(CO)C)C)O)C)[C@H](C2)C '
               'NAME: Neomangicol A REASON: MISSED Missing required pyran or '
               'modified oxygen heterocycle\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Insufficient oxygen atoms for '
                                               'iridoid structure'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Molecular weight outside '
                                               'typical range for iridoids'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Molecular weight outside '
                                               'typical range for iridoids'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Insufficient oxygen atoms for '
                                               'iridoid structure'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Missing characteristic '
                                               'iridoid/secoiridoid skeleton'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Missing characteristic '
                                               'iridoid/secoiridoid skeleton'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Insufficient oxygen atoms for '
                                               'iridoid structure'},
                                 {   'smiles': 'C1CN=C(Cc2ccccc2)N1',
                                     'name': 'tolazoline',
                                     'reason': 'Insufficient oxygen atoms for '
                                               'iridoid structure'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)C3=CC=NC=C3)O[C@@H]1CN(C)S(=O)(=O)C)[C@@H](C)CO',
                                     'name': 'N-[(2S,3R)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl(methylsulfonyl)amino]methyl]-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-8-yl]-4-pyridinecarboxamide',
                                     'reason': 'Missing characteristic '
                                               'iridoid/secoiridoid skeleton'},
                                 {   'smiles': 'CC(C)[C@@H]1CC\\C(C)=C\\CC[C@@](C)(O)\\C=C\\1',
                                     'name': 'germacra-1(10),5-dien-4alpha-ol',
                                     'reason': 'Insufficient oxygen atoms for '
                                               'iridoid structure'}],
    'sample_false_negatives': [   {   'smiles': 'OC12C(C(CC1(CO)C)(C)C)CCC(C2C(O)=O)C',
                                      'name': '3a-hydroxy-3-(hydroxymethyl)-1,1,3,5-tetramethyl-2,4,5,6,7,7a-hexahydroindene-4-carboxylic '
                                              'acid',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'C12=C3[C@H]([C@@H]4C(CC[C@@]4(C3)C)(C)C)CC[C@]1([C@H](C)CC2)C',
                                      'name': 'Cattleyene',
                                      'reason': 'Insufficient oxygen atoms for '
                                                'iridoid structure'},
                                  {   'smiles': 'C[C@@H]1C[C@@]2(CC(C)=CCC[C@@H]2[C@@]1(CO)CC\\C(CO)=C\\CO)C=O',
                                      'name': 'Portulal',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'O=C1NCCC[C@@H]2NC(=O)C(C2=O)=C(O)C=C[C@H]3[C@H](CC=C1)[C@H]4[C@H](O)[C@@H](O)[C@H]5[C@H]([C@H](O)C)[C@@H](C[C@@H]5[C@H]4C3)C',
                                      'name': 'Hydroxyikarugamycin C',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'C1(=C[C@H]2[C@@]3([C@H](CC[C@@H]3C2(C)C)C)CC1)C',
                                      'name': 'Italicene',
                                      'reason': 'Insufficient oxygen atoms for '
                                                'iridoid structure'},
                                  {   'smiles': 'Cl/C=C\\1/C2=C3C(=C[C@]4(C[C@@H]([C@@]([C@H]4[C@@H]3CC1)(C[C@H](O)[C@H](O)[C@@](O)(CO)C)C)O)C)[C@H](C2)C',
                                      'name': 'Neomangicol A',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'O=C1N[C@@H](CCC(=O)O)C(C1=C(O)/C=C/[C@H]2[C@@H]([C@H]3C=C[C@@H]4[C@@H]([C@H]3C2)C[C@@H](C)[C@H]4C(=O)C)C/C=C\\C(=O)N)=O',
                                      'name': 'Clifednamide D',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'O=C1N[C@@H](CCC(=O)O)C(C1=C(O)/C=C/C=C/[C@H]2[C@H](C=C[C@@H]3[C@@H]2C[C@@H](C)[C@H]3CC)/C=C/C=C/C(=O)N)=O',
                                      'name': 'Clifednamide H',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'O[C@H]1[C@@H](C(=C(C1)C)CO)CCO',
                                      'name': '1-Deoxyeucommiol',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid skeleton'},
                                  {   'smiles': 'O=C1C=C2C3=C(O[C@]4([C@@]2([C@@H](C(=C)C)CC4)CC1)C)C=CC(=C3)O',
                                      'name': 'Ganocin B',
                                      'reason': 'Missing characteristic '
                                                'iridoid/secoiridoid '
                                                'skeleton'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code based on '
                 'the outcomes:\n'
                 '\n'
                 '1. False Negatives Analysis:\n'
                 '- Many legitimate iridoids were missed because the '
                 'pyran/heterocycle detection was too strict\n'
                 '- The code failed to identify modified or rearranged iridoid '
                 'skeletons (like in Genipin)\n'
                 '- Some compounds with fused ring systems were missed\n'
                 '- The cyclopentane/seco pattern matching was too '
                 'restrictive\n'
                 '\n'
                 '2. False Positives Analysis:\n'
                 '- Many glucosides and other sugar-containing compounds were '
                 'wrongly classified because they contain pyran rings and '
                 'similar substituents\n'
                 '- Complex molecules with multiple rings and '
                 'oxygen-containing heterocycles were misidentified\n'
                 '- The presence of hydroxyl, carboxyl, and methyl groups '
                 "alone wasn't specific enough\n"
                 '\n'
                 "Here's how we can improve the code:\n"
                 '\n'
                 '1. Make the core structure detection more flexible:\n'
                 '- Use multiple SMARTS patterns for different iridoid '
                 'skeleton variations\n'
                 '- Include patterns for rearranged and modified iridoid '
                 'skeletons\n'
                 '- Make the cyclopentane detection more flexible to catch '
                 'fused systems\n'
                 '\n'
                 '2. Add more specific criteria:\n'
                 '- Check for carbon skeleton connectivity patterns specific '
                 'to iridoids\n'
                 '- Add molecular weight range check\n'
                 '- Implement better discrimination between iridoid-specific '
                 'oxygen heterocycles and sugar rings\n'
                 '- Add checks for specific carbon frameworks common in '
                 'iridoids\n'
                 '\n'
                 '3. Improve substituent analysis:\n'
                 '- Consider relative positions of substituents\n'
                 '- Check for specific substitution patterns characteristic of '
                 'iridoids\n'
                 '- Add more specific patterns for iridoid-specific functional '
                 'group arrangements\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 72,
    'num_false_positives': 100,
    'num_true_negatives': 614,
    'num_false_negatives': 37,
    'num_negatives': None,
    'precision': 0.4186046511627907,
    'recall': 0.6605504587155964,
    'f1': 0.5124555160142349,
    'accuracy': 0.8335358444714459,
    'negative_predictive_value': 0.9431643625192012}