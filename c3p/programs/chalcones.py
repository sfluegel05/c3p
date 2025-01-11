"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone),
    ArC(=O)CH=CHAr, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chalcone core pattern: aromatic ring connected via α,β-unsaturated ketone to another aromatic ring
    chalcone_pattern = Chem.MolFromSmarts("[a][C](=O)[C]=[C][a]")
    # Define dihydrochalcone pattern: saturated version of chalcone
    dihydrochalcone_pattern = Chem.MolFromSmarts("[a][C](=O)[C][C][a]")

    # Check for the chalcone or dihydrochalcone core
    is_chalcone = mol.HasSubstructMatch(chalcone_pattern)
    is_dihydrochalcone = mol.HasSubstructMatch(dihydrochalcone_pattern)

    if not (is_chalcone or is_dihydrochalcone):
        return False, "Chalcone core not found"

    # Verify the number of aromatic rings (should be at least 2)
    aromatic_rings = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            aromatic_rings.append(ring)
    num_aromatic_rings = len(aromatic_rings)
    if num_aromatic_rings < 2:
        return False, f"Less than two aromatic rings detected ({num_aromatic_rings})"

    # Find all matches to the chalcone or dihydrochalcone pattern
    if is_chalcone:
        matches = mol.GetSubstructMatches(chalcone_pattern)
        pattern_name = "chalcone"
    else:
        matches = mol.GetSubstructMatches(dihydrochalcone_pattern)
        pattern_name = "dihydrochalcone"

    if matches:
        return True, f"{pattern_name.capitalize()} core detected"
    else:
        return False, "Chalcone core not found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23086',
                          'name': 'chalcones',
                          'definition': 'A ketone that is '
                                        '1,3-diphenylpropenone '
                                        '(benzylideneacetophenone), '
                                        'ArCH=CH(=O)Ar, and its derivatives '
                                        'formed by substitution.',
                          'parents': [   'CHEBI:51689',
                                         'CHEBI:72544',
                                         'CHEBI:76224'],
                          'xrefs': ['KEGG:C15589'],
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
               'COc1ccc(\\C=C\\C(=O)c2cc(CC=C(C)C)c(OC)cc2O)cc1 NAME: '
               "4'-O-methylbavachalcone REASON: MISSED Chalcone core not found "
               'or is part of a ring system\n'
               ' * SMILES: B(C1=CC=C(C=C1)C(=O)C=CC2=CC=C(C=C2)I)(O)O NAME: '
               '[4-[3-(4-iodophenyl)prop-2-enoyl]phenyl]boronic acid REASON: '
               'MISSED Chalcone core not found or is part of a ring system\n'
               ' * SMILES: COc1c(C)c(O)c(C)c(O)c1C(=O)\\C=C\\c1ccc(O)cc1 NAME: '
               "(E)-4,2',4'-trihydroxy-6'-methoxy-3',5'-dimethylchalcone "
               'REASON: MISSED Chalcone core not found or is part of a ring '
               'system\n'
               ' * SMILES: O(C1=C(/C=C/C(=O)C2=CC=C(OC)C=C2)C=CC(O)=C1)C NAME: '
               'Glypallichalcone REASON: MISSED Chalcone core not found or is '
               'part of a ring system\n'
               ' * SMILES: O(C=1C(=C(O)C=C(O)C1)C(=O)/C=C/C2=CC=CC=C2)C NAME: '
               'Cardamonin REASON: MISSED Chalcone core not found or is part '
               'of a ring system\n'
               ' * SMILES: '
               'O1C(C=CC=2C1=C(C(O)=CC2OC)C(=O)/C=C/C3=CC=C(O)C=C3)(C)C NAME: '
               'Citrunobin REASON: MISSED Chalcone core not found or is part '
               'of a ring system\n'
               ' * SMILES: '
               'O1C2=C(C(=C1O)C(=O)/C=C/C3=CC=CC=C3)C(OC)=C(OC)C(O)=C2C(=O)/C=C/C4=CC=CC=C4 '
               'NAME: Didymocalyxin B REASON: MISSED Chalcone core not found '
               'or is part of a ring system\n'
               ' * SMILES: COC1=CC(=CC(=C1OC)OC)C=CC(=O)C2=CC=CC=C2 NAME: '
               '1-phenyl-3-(3,4,5-trimethoxyphenyl)-2-propen-1-one REASON: '
               'MISSED Chalcone core not found or is part of a ring system\n'
               ' * SMILES: O1C(C=CC2=C1C=CC(=C2O)C(=O)/C=C/C3=CC=CC=C3)(C)C '
               'NAME: Lonchocarpin REASON: MISSED Chalcone core not found or '
               'is part of a ring system\n'
               ' * SMILES: C1=CC(=CC=C1C=CC(=O)C2=C(C=C(C=C2)O)O)O NAME: '
               '1-(2,4-dihydroxyphenyl)-3-(4-hydroxyphenyl)-2-propen-1-one '
               'REASON: MISSED Chalcone core not found or is part of a ring '
               'system\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Chalcone core not found'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Chalcone core not found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'The previous code failed to correctly identify chalcones '
                 'because the SMARTS pattern used to define the chalcone core '
                 'was incorrect. Specifically, the pattern '
                 '`[a]-[C;R0](=O)-[C;R0]=[C;R0]-[a]` does not accurately '
                 'represent the connectivity of atoms in a chalcone.\n'
                 '\n'
                 'A chalcone consists of two aromatic rings connected via an '
                 'α,β-unsaturated ketone. The correct core structure is:\n'
                 '\n'
                 '- An aromatic ring connected to a carbonyl group (C=O)\n'
                 '- The carbonyl carbon is connected to an alkene (C=C)\n'
                 '- The alkene is connected to another aromatic ring\n'
                 '\n'
                 'This can be represented by the SMARTS pattern '
                 '`[a][C](=O)[C]=[C][a]`. Additionally, excluding chalcones '
                 'where the core is part of a ring system is too restrictive, '
                 'as many chalcone derivatives have substitutions that form '
                 'rings involving the core atoms.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '\n'
                 '1. Correct the SMARTS pattern to accurately represent the '
                 'chalcone core.\n'
                 '2. Remove the restriction that excludes molecules where the '
                 'chalcone core is part of a ring system.\n'
                 '\n'
                 'Here is the updated code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 223,
    'num_false_positives': 100,
    'num_true_negatives': 12036,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6904024767801857,
    'recall': 1.0,
    'f1': 0.8168498168498168,
    'accuracy': 0.9919087304798123,
    'negative_predictive_value': 1.0}