"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenes, containing oxygen atoms in addition to the typical carotenoid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for the presence of a long polyene chain (carotenoid backbone)
    # Use a more flexible pattern to identify the carotenoid backbone
    carotenoid_pattern = Chem.MolFromSmarts("[CH2,CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
        return False, "No carotenoid backbone found"

    # Check for specific oxygen-containing functional groups (hydroxyl, epoxy, carbonyl, ether, carboxyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("[OX2]C1CC1")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    ether_pattern = Chem.MolFromSmarts("[OX2]C")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_hydroxyl or has_epoxy or has_carbonyl or has_ether or has_carboxyl):
        return False, "No hydroxyl, epoxy, carbonyl, ether, or carboxyl groups found"

    # Additional check for conjugated double bonds
    conjugated_double_bonds = Chem.MolFromSmarts("[CH,CH0]=[CH,CH0]~[CH,CH0]=[CH,CH0]")
    if not mol.HasSubstructMatch(conjugated_double_bonds):
        return False, "No conjugated double bonds found, not a carotenoid"

    return True, "Contains a carotenoid backbone with oxygenated functional groups and conjugated double bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27325',
                          'name': 'xanthophyll',
                          'definition': 'A subclass of carotenoids consisting '
                                        'of the oxygenated carotenes.',
                          'parents': ['CHEBI:23044'],
                          'xrefs': ['DrugBank:DB00137'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7682119205298014 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C(O[C@H](C(O)(C)C)CC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C[C@@H](O)CC2(C)C)C)\\C)\\C)/C)/C)/C)[C@@H](O)[C@H](O)[C@@H]([C@@H]1C)O '
               "NAME: (3R,2'S)-myxol 2'-fucoside REASON: CORRECT Contains a "
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O[C@H]1CC([C@](C1)(C)C(=O)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=2C(CCCC2C)(C)C)/C)/C)(C)C '
               'NAME: Cryptocapsin REASON: CORRECT Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\[C@H]1C(C)=C[C@H](O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C '
               'NAME: alpha-cryptoxanthin REASON: CORRECT Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'COC(C)(C)C(=O)\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CC\\C=C(/C)CCCC(C)(C)O '
               'NAME: Hydroxyspheroidenone REASON: CORRECT Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'OC1CC(C(C(=C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C2C(CC(=O)C=C2C)(C)C)/C)/C)(C)C '
               "NAME: 3'-hydroxy-e,e-caroten-3-one REASON: CORRECT Contains a "
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O1C2(C1(CCCC2(C)C)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(CC(O)CC3C)(C)C)\\C)\\C '
               "NAME: 5',6'-Epoxy-5',6'-dihydro-beta,beta-caroten-3-ol REASON: "
               'CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C(O)C(O[C@@H]2OC[C@@H](O)[C@@H]([C@H]2O)O)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)CC1O)C '
               "NAME: 1'-beta-xylosyl-3,2'-dihydroxy-beta,psi-caroten-4-one "
               'REASON: CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'S(O[C@H]1CC(C(=C(C1=O)C)/C=C/C(/C)=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=2C([C@@H](O)[C@H](O)CC2C)(C)C)\\C)\\C)/C)(C)C)(O)(=O)=O '
               'NAME: Erythroxanthin sulfate REASON: CORRECT Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C1=C(C)C(=O)C(=O)CC1(C)C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C(=O)C(=O)CC1(C)C '
               'NAME: astacene REASON: CORRECT Contains a carotenoid backbone '
               'with oxygenated functional groups\n'
               ' * SMILES: '
               'O1C2(C1(CCCC2(C)C)C)\\C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C/C=C(/C=C/C=C(\\C3OC4(C(C(CCC4)(C)C)=C3)C)/C)\\C '
               "NAME: (5R,5'R,6S,8'R)-Luteochrome REASON: CORRECT Contains a "
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'C1(=C(C[C@H](CC1(C)C)O)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/[C@H]2C(C=CCC2(C)C)=C)\\C)\\C)/C)/C '
               "NAME: (3R,6'R)-3-hydroxy-3',4'-didehydro-beta,gamma-carotene "
               'REASON: CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C\\CC(O)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1O)C '
               'NAME: Flexixanthin REASON: CORRECT Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C(OC[C@H]1O[C@@H](OC(C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/CC(O)(C)C)\\C)\\C)\\C)/C)/C)/C)(C)C)[C@H](O)[C@H]([C@@H]1O)O)CCCCCCCCCCCCC '
               'NAME: I\' -(6"-0-myristyl-p-glucopyanosyloxy) -3,4,3\', 4\' '
               "-tetradehydro1,2,1',2'-tetrahydro-l-hydroxy-yl,yl-carotene "
               'REASON: CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C2=C(C)C[C@@H](O)CC2(C)C)C(C)(C)O)[C@@H](O)[C@H](O)[C@@H]1O '
               "NAME: (3R,2'S)-Myxol 2'-alpha-L-fucoside REASON: CORRECT "
               'Contains a carotenoid backbone with oxygenated functional '
               'groups\n'
               ' * SMILES: '
               'S(O[C@H]1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C([C@@H](O)[C@H](O)CC2C)(C)C)\\C)\\C)(C)C)(O)(=O)=O '
               'NAME: Caloxanthin sulfate REASON: CORRECT Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'OC1CC(C(=C([C@H]1O)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=2C(CC(O)C(=O)C2C)(C)C)\\C)\\C)(C)C '
               'NAME: Idoxanthin REASON: CORRECT Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1CCC(C(=C1C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/[C@H](O)C(O)(C)C)\\C)\\C)\\C)(C)C '
               "NAME: 1',2'-Dihydro-1',2'-dihydroxy-4-ketotoruelene/ "
               '(4-Ketophleixanthophyll) REASON: CORRECT Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CC(\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\[C@H](O)C(C)(C)O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)C(=O)[C@@H](O)CC1(C)C '
               'NAME: 4-Ketomyxol REASON: CORRECT Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O([C@H]1CC(C(=C(C1)C)\\C=C\\C(\\C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=2C(C[C@H](O)C(=O)C2C)(C)C)/C)/C)(C)C)[C@@H]3O[C@@H]([C@@H](O)C(O)C3O)CO '
               "NAME: Adonixanthin 3'-beta-D-glucoside/ Adonixanthin "
               "3'-glucoside REASON: CORRECT Contains a carotenoid backbone "
               'with oxygenated functional groups\n'
               ' * SMILES: '
               'O[C@@H]1CC([C@@H](C(=C1)C)/C=C/C(/C)=C\\C=C\\C(\\C)=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C=2C(C[C@H](O)CC2C)(C)C)/C)/C)(C)C '
               "NAME: (3R,3'R,6'R,9'-cis)-b,e-Carotene-3,3'-diol REASON: "
               'CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'C\\C(CCCC(C)(C)O)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(\\C)CCCC(C)(C)O '
               "NAME: 1,1'-dihydroxy-1,1',2,2'-tetrahydrolycopene REASON: "
               'CORRECT Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/[C@]23O[C@@]2(C)CCCC3(C)C)C(C)(C)CCC1 '
               'NAME: (5S,6R)-beta-carotene 5,6-epoxide REASON: CORRECT '
               'Contains a carotenoid backbone with oxygenated functional '
               'groups\n'
               ' * SMILES: '
               '[H]C(=C=C1C(C)(C)C[C@H](O)C[C@@]1(C)O)C(\\C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)/C=C/[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C '
               "NAME: 9'-cis-neoxanthin REASON: CORRECT Contains a carotenoid "
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1C(=C(/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/[C@H](OC)C(O)(C)C)\\C)\\C)\\C)/C)/C)C(C)(C)C[C@@H]1O)C '
               "NAME: 2'-Methoxyflexixanthin REASON: CORRECT Contains a "
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'C(/C=C/C(=C/C=C/C(=C/C=C/C=C(\\C)/C=C/C=C(/C=C/[C@@]12[C@](C[C@@H](O)CC1(C)C)(C)O2)\\C)/C)/C)([C@@]3(C[C@H](CC3(C)C)O)C)=O '
               'NAME: capsanthin 5,6-epoxide REASON: CORRECT Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               'False positives: SMILES: '
               'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C '
               'NAME: Concanamycin D REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C '
               'NAME: 2,2-dimethylbutanoic acid '
               '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
               'ester REASON: WRONGLY CLASSIFIED Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)[C@H]5OC4)COC(=O)C=6NC=CC6)O)C)C)O[C@@H]([C@@H](C)[C@@H]([C@H]3O)OC(=O)C(CC(C)C)C)C '
               'NAME: Milbemycin alpha29 REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'ClC1=C(OC)C(=C(C)C(=C1O)Cl)C(=O)O[C@@]2(C(=O)C=3C(C=C(/C=C/C)OC3)=CC2=O)C '
               'NAME: Falconensin H REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@]2(C[C@H](C1)O)[H])[H])(CC[C@@H]4O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)C(=O)[O-])O)O)O)[H])C)[H])C '
               'NAME: 5alpha-androstane-3beta,17beta-diol '
               '17-O-(beta-D-glucuronide)(1-) REASON: WRONGLY CLASSIFIED '
               'Contains a carotenoid backbone with oxygenated functional '
               'groups\n'
               ' * SMILES: '
               '[C@]12(C=CC3=C(C=C(C(=C3C)O)O)C1=CC[C@@]4([C@@]2(CC[C@@]5([C@]4(C[C@@](CC5)(C(O)=O)C)[H])C)C)C)C '
               'NAME: wilforol B REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O(C1C(C2C(C3C(C4(C(CC3)=C5C(CC4)(CCC(C5)(C)C)C(OC)=O)C)(CC2)C)(CC1)C)(C)C)C(=O)C '
               'NAME: Methyl 3b-hydroxy-13(18)-oleanen-28-oate REASON: WRONGLY '
               'CLASSIFIED Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'O=C1OC(C2=CN=CC=C2)=CC3=C1C(O)C4C5(C(C(C(OC(=O)C)CC5)(COC(=O)C)C)CC(C4(O3)C)OC(=O)CC)C '
               'NAME: Pyripyropene C REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C(O)[C@@H]([C@@H]1[C@@]2([C@@](C=3C([C@@]4(C(C([C@@H](O)CC4)(C)C)CC3)C)=CC2)(C)CC1)C)CCC=C(C)C '
               'NAME: 3beta-Hydroxylanosta-7, 9 (11), 24-trien-21-oic acid '
               'REASON: WRONGLY CLASSIFIED Contains a carotenoid backbone with '
               'oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1C=2OC([C@H]3C[C@H](OC(=O)C)C(=C[C@H]3C2C(=O)C=4C1=C(O)C(=C(O)C4)C)C)(C)C '
               'NAME: Unnamed naphterpin 2 REASON: WRONGLY CLASSIFIED Contains '
               'a carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O1C[C@@]2([C@]3([C@@]([C@@]4([C@@](CC3)(C([C@@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)CC4)(C)C)[H])C)(CC[C@@]2([C@H]([C@@H](OC(=O)C)C[C@H]6OC(=O)[C@@H](C6)C)C)[H])[H])C)CC1=O '
               'NAME: Hovenidulcioside B2 REASON: WRONGLY CLASSIFIED Contains '
               'a carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O(C1CC=2C(C3C(C4C(C(CC4)C(CCC(C(C)C)C)C)(CC3)C)CC2)(CC1)C)C(=O)CCCCCCC/C=C\\CCCCCCCC '
               'NAME: Campesteryl elaidate REASON: WRONGLY CLASSIFIED Contains '
               'a carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O1C2C(C3(C(C4C(C5(C(CC4)CC(OC6OC(C(O)C(OC7OC(C(O)C(O)C7O)C)C6OC8OC(C(O)C(O)C8O)CO)CO)CC5)C)CC3)C2)C)C(C1(O)CCC(COC9OC(C(O)C(O)C9O)CO)C)C '
               'NAME: Melongoside P REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@]3([H])[C@]([H])(CC(=O)[C@@]4(C)[C@@]3([H])C[C@]3([H])O[C@]5(CC[C@@H](C)CO5)[C@@H](C)[C@]43[H])[C@@]1(C)CC[C@@H](C2)O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O '
               'NAME: agavoside A REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'C(CCl)(=O)[C@]12[C@@]3([C@@](C[C@]1(OC(O2)(C)C)[H])([C@]4([C@@]([C@H](C3)O)(F)[C@]5(C)C(CC4)=CC(CC5)=O)[H])[H])C '
               'NAME: Halcinonide REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CC1=CC(=O)[C@@H](O)[C@@]2(C)[C@H]1C[C@H]1OC(=O)[C@H](O)[C@@]3(O)[C@]4(C)OC[C@]13[C@@H]2[C@@H](O)[C@@H]4O '
               'NAME: bruceine D REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C(O[C@@H]1OC[C@@H](O)[C@@H]([C@H]1OC(=O)CCCCCCCCCCC(C)C)O)/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C(=O)O)\\C)\\C)\\C)/C)/C)/C '
               'NAME: Diapolycopenedioic acid xylosyl ester REASON: WRONGLY '
               'CLASSIFIED Contains a carotenoid backbone with oxygenated '
               'functional groups\n'
               ' * SMILES: '
               'O=C1C2=C([C@@]3([C@@H](O)C[C@]4([C@]3(C1)C)O[C@]5(OC(=O)[C@H](C5)C)C[C@H]4C)C)[C@@H](OC(=O)C)C[C@@H]6[C@@]2(CCC(C6(C)C)=O)C '
               'NAME: Spiroganocalitone D REASON: WRONGLY CLASSIFIED Contains '
               'a carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C=CC=C2C3=C4C=5C(C(=O)C=C4)=C(O)C=C(C5[C@@]63[C@H](C1)CC(=O)C7=C6C=C(O)C=C7O)O '
               'NAME: Dalesconol C REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@@](C)(OC(C)=O)[C@@H]2O[C@H](C[C@]1(C)O)[C@@H]1[C@H]2[C@@H](C(C)C)[C@H](OC(C)=O)[C@H](OC(C)=O)[C@@]11CO1 '
               'NAME: klymollin G REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C(O[C@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1O)O)O)[C@@H]([C@@H]2[C@@]3([C@@](C4=C([C@@]5([C@H](C([C@H](O)CC5)(C)C)CC4)C)CC3)(C)CC2)C)CCC=C(C)C '
               'NAME: Fomitoside I REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O(C1C2C(C3C(C1O)C(C=4C(OC(=O)C4)(O)C3)C)(CCCC2(CO)C)C)C(=O)C '
               'NAME: '
               '[6,10a-dihydroxy-4-(hydroxymethyl)-4,7,11b-trimethyl-9-oxo-1,2,3,4a,5,6,6a,7,11,11a-decahydronaphtho[2,1-f][1]benzofuran-5-yl] '
               'acetate REASON: WRONGLY CLASSIFIED Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O=C1C2=C([C@]3(CCC(C([C@@H]3C1)(C)C)=O)C)CC[C@]4([C@]2(CC[C@@H]4[C@@H](C[C@H]5OC(=O)C(=C5)C)C)C)C '
               'NAME: '
               '(+)-(5alpha,23R,24Z)-lanosta-8,24-dien-3,7-dioxo-23,26-gamma-lactone '
               'REASON: WRONGLY CLASSIFIED Contains a carotenoid backbone with '
               'oxygenated functional groups\n'
               ' * SMILES: '
               'CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O[C@H](C)CNC(=O)CC[C@@]4([C@H]([C@@H]5[C@]6([C@@]([C@@H](C(=C(C7=NC(=CC8=NC(=C(C4=N5)C)[C@H](C8(C)C)CCC(=O)N)[C@H]([C@]7(C)CC(=O)N)CCC(=O)N)C)[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3] '
               'NAME: cobalt(3+) '
               '[(2R,3S,4R,5S)-5-(5,6-dimethylbenzimidazol-1-yl)-4-hydroxy-2-(hydroxymethyl)oxolan-3-yl] '
               '[(2R)-1-[3-[(1R,2R,3R,4Z,7S,9Z,12S,13S,14Z,17S,18S,19R)-2,13,18-tris(2-amino-2-oxoethyl)-7,12,17-tris(3-amino-3-oxopropyl)-3,5,8,8,13,15,18,19-octamethyl-2,7,12,17-tetrahydro-1H-corrin-21-id-3-yl]propanoylamino]propan-2-yl] '
               'hydrogen phosphate cyanidecobalt(3+) '
               '[(2R,3S,4R,5S)-5-(5,6-dimethylbenzimidazol-1-yl)-4-hydroxy-2-(hydroxymethyl)oxolan-3-yl] '
               '[(2R)-1-[3-[(1R,2R,3R,5Z,7S,10Z,12S,13S,15Z,17S,18S,19R)-2,13,18-tris(2-amino-2-oxoethyl)-7,12,17-tris(3-amino-3-oxopropyl)-3,5,8,8,13,15,18,19-octamethyl-2,7,12,17-tetrahydro-1H-corrin-24-id-3-yl]propanoylamino]propan-2-yl] '
               'hydrogen phosphate cyanidecobalt(3+) '
               '[(2R,3S,4R,5S)-5-(5,6-dimethylbenzimidazol-1-yl)-4-hydroxy-2-(hydroxymethyl)oxolan-3-yl] '
               '[(2R)-1-[3-[(1R,2R,3R,7S,12S,13S,17S,18S,19R)-2,13,18-tris(2-amino-2-oxoethyl)-7,12,17-tris(3-amino-3-oxopropyl)-3,5,8,8,13,15,18,19-octamethyl-2,7,12,17-tetrahydro-1H-corrin-24-id-3-yl]propanoylamino]propan-2-yl] '
               'hydrogen REASON: WRONGLY CLASSIFIED Contains a carotenoid '
               'backbone with oxygenated functional groups\n'
               ' * SMILES: '
               'O1C([C@@H]2[C@@]3(O)C=4C(=CC=C5C4C=6[C@@H]1O[C@]7(OC)[C@]([C@@]8(C([C@@]9%10O[C@@H]9[C@@H](O)[C@H](OC%10CC8)C(=C)C)CC7)C)(C)C6N5)CC([C@H]3C2)=C)(C)C '
               'NAME: Rhizovarin C REASON: WRONGLY CLASSIFIED Contains a '
               'carotenoid backbone with oxygenated functional groups\n'
               'False negatives: SMILES: '
               'CC(\\C=C\\C=C(/C)C(O)=O)=C/C=C/C=C(C)/C=C/C=C(\\C)C(O)=O NAME: '
               'crocetin REASON: MISSED Molecular weight too low for a '
               'xanthophyll\n'
               ' * SMILES: '
               'C\\C(C=O)=C/C=C/C=C(\\C)/C=C/C=C(\\C)C=[C@]=C1C(C)(C)C[C@H](O)C[C@@]1(C)O '
               'NAME: '
               "(3S,5R,6R)-3,5-dihydroxy-6,7-didehydro-5,6-dihydro-12'-apo-beta-caroten-12'-al "
               'REASON: MISSED Molecular weight too low for a xanthophyll\n'
               ' * SMILES: '
               'C\\C(C=O)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C '
               'NAME: '
               "(3S,5R,6S)-5,6-epoxy-3-hydroxy-5,6-dihydro-12'-apo-beta-caroten-12'-al "
               'REASON: MISSED Molecular weight too low for a xanthophyll\n'
               ' * SMILES: '
               'C(\\C=C\\C(OC)=O)(=C\\C=C\\C(=C\\C=C\\C=C(\\C=C\\C=C(\\C=C\\C(=O)O)/C)/C)\\C)/C '
               'NAME: bixin REASON: MISSED Molecular weight too low for a '
               'xanthophyll\n'
               ' * SMILES: '
               '[H]C(=O)\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C1=C(C)CCCC1(C)C '
               "NAME: 14'-apo-beta-carotenal REASON: MISSED Molecular weight "
               'too low for a xanthophyll\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No carotenoid backbone found'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'No oxygen atoms found, not a '
                                               'xanthophyll'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No conjugated double bonds '
                                               'found, not a carotenoid'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'The previous program had a few issues that led to incorrect '
                 'classifications:\n'
                 '\n'
                 '1. **Molecular Weight Threshold**: The molecular weight '
                 'threshold of 400 Da was too restrictive, causing some valid '
                 'xanthophylls (like crocetin) to be incorrectly classified as '
                 'non-xanthophylls. Xanthophylls can have a wide range of '
                 'molecular weights, and some smaller molecules (like '
                 'crocetin) are still valid xanthophylls.\n'
                 '\n'
                 '2. **Carotenoid Backbone Pattern**: The pattern used to '
                 'identify the carotenoid backbone was too rigid and might not '
                 'capture all variations of the carotenoid structure. The '
                 'pattern should be more flexible to account for different '
                 'degrees of unsaturation and branching.\n'
                 '\n'
                 '3. **Oxygen-Containing Functional Groups**: The program '
                 'checked for specific oxygen-containing functional groups, '
                 'but the patterns used were not comprehensive enough. For '
                 'example, the epoxy pattern (`[OX2]C1CC1`) might miss some '
                 'epoxy groups that are part of larger ring systems.\n'
                 '\n'
                 '4. **False Positives**: Some molecules that are not '
                 'xanthophylls were incorrectly classified as xanthophylls '
                 'because they contained oxygen atoms and had a molecular '
                 'weight above 400 Da. The program needs to better distinguish '
                 'between xanthophylls and other oxygen-containing compounds.\n'
                 '\n'
                 '### Improvements:\n'
                 '1. **Relax the Molecular Weight Threshold**: Remove or lower '
                 'the molecular weight threshold to allow smaller xanthophylls '
                 'to be correctly classified.\n'
                 '\n'
                 '2. **Refine the Carotenoid Backbone Pattern**: Use a more '
                 'flexible pattern to capture the carotenoid backbone, '
                 'considering different degrees of unsaturation and '
                 'branching.\n'
                 '\n'
                 '3. **Improve Oxygen-Containing Functional Group Detection**: '
                 'Use more comprehensive patterns to detect oxygen-containing '
                 'functional groups, including those in more complex ring '
                 'systems.\n'
                 '\n'
                 '4. **Add Additional Checks**: Include checks for specific '
                 'structural features that are common in xanthophylls, such as '
                 'the presence of conjugated double bonds and specific '
                 'oxygen-containing groups.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 179,
    'num_false_positives': 100,
    'num_true_negatives': 3154,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6415770609318996,
    'recall': 1.0,
    'f1': 0.7816593886462881,
    'accuracy': 0.9708709583454704,
    'negative_predictive_value': 1.0}