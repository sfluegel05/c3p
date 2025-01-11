"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:8338 porphyrins
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a fundamental skeleton of four pyrrole nuclei united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate if porphyrin ring is found
    porphyrin_found = False

    # Iterate over all rings
    for ring in atom_rings:
        # Check if ring is 16-membered
        if len(ring) == 16:
            # Count the number of nitrogen atoms in the ring
            num_nitrogens = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'N':
                    num_nitrogens += 1
            # Check if there are exactly 4 nitrogen atoms
            if num_nitrogens == 4:
                porphyrin_found = True
                break  # No need to check other rings

    if porphyrin_found:
        return True, "Porphyrin ring system detected (16-membered ring with 4 nitrogen atoms)"
    else:
        return False, "No porphyrin ring system detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26214',
                          'name': 'porphyrins',
                          'definition': 'Natural pigments containing a '
                                        'fundamental skeleton of four pyrrole '
                                        'nuclei united through the '
                                        'alpha-positions by four methine '
                                        'groups to form a macrocyclic '
                                        'structure.',
                          'parents': ['CHEBI:36309'],
                          'xrefs': ['Wikipedia:Porphyrin'],
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
               'CC[C@@H]1[C@@H](C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(O)=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31 '
               'NAME: 3-vinylbacteriochlorophyllide a REASON: MISSED No '
               'porphyrin ring system detected\n'
               ' * SMILES: '
               'COC(=O)[C@H]1C(=O)c2c(C)c3C=C4C(C=C)=C(C)C5=[N+]4[Mg--]46n3c2C1=C1[C@@H](CCC(O)=O)[C@H](C)C(C=c2c(C)c(C=C)c(=C5)n42)=[N+]61 '
               'NAME: divinyl chlorophyllide a REASON: MISSED No porphyrin '
               'ring system detected\n'
               ' * SMILES: '
               'COC(=O)CCc1c(C)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CCC(O)=O)c5C)c(C=C)c4C)C1=CC=C([C@H](C(=O)OC)[C@]31C)C(=O)OC '
               'NAME: '
               '(2S,2(1)R)-8-ethenyl-2(1),2(2)-bis(methoxycarbonyl)-17-(3-methoxy-3-oxopropyl)-2,7,12,18-tetramethyl-2,2(1)-dihydrobenzo[b]porphyrin-13-propanoic '
               'acid REASON: MISSED No porphyrin ring system detected\n'
               ' * SMILES: '
               'OC(=O)CCc1c(CC(O)=O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CCC(O)=O)c5CC(O)=O)c(CCC(O)=O)c4CC(O)=O)c(CCC(O)=O)c3CC(O)=O '
               'NAME: uroporphyrin III REASON: MISSED No porphyrin ring system '
               'detected\n'
               ' * SMILES: Bc1c2ccc(cc3ccc(cc4ccc(cc5ccc1[nH]5)n4)[nH]3)n2 '
               'NAME: 5-borylporphyrin REASON: MISSED No porphyrin ring system '
               'detected\n'
               ' * SMILES: '
               'N1C2=C(C=3NC(C(=C4N=C(C(=C5N=C(C(=C1C=C2)C6=CC=NC=C6)C=C5)C7=CC=NC=C7)C=C4)C8=CC=NC=C8)=CC3)C9=CC=NC=C9 '
               'NAME: 5,10,15,20-tetra(4-pyridyl)-21H,23H-porphine REASON: '
               'MISSED No porphyrin ring system detected\n'
               ' * SMILES: '
               'Cc1c(C=C)c2C=C3[N+]4=C(C=c5c(CCC(O)=O)c(C)c6=CC7=[N+]8C(=Cc1n2[Fe--]48n56)C(C=C)=C7C)[C@](O)(CCC(O)=O)[C@@]3(C)O '
               'NAME: heme d cis-diol REASON: MISSED No porphyrin ring system '
               'detected\n'
               ' * SMILES: '
               'CC1OC(=O)C[C@@]2(C)[C@H](CCC(O)=O)C3=CC4=[N+]5C(=CC6=[N+]7C(=CC8=[N+]9C(=C(CC(O)=O)[C@@]8(C)CCC(O)=O)C12N3[Co--]579)C(CCC(O)=O)=C6CC(O)=O)[C@@H](CCC(O)=O)[C@]4(C)CC(O)=O '
               'NAME: cobalt(II)-factor IV REASON: MISSED No porphyrin ring '
               'system detected\n'
               ' * SMILES: '
               'Cc1c(CCC(O)=O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CCC(O)=O)c5C)c(CCC(O)=O)c4C)c(CCC(O)=O)c3C '
               'NAME: coproporphyrin I REASON: MISSED No porphyrin ring system '
               'detected\n'
               ' * SMILES: O=Cc1cc2cc3ccc(cc4ccc(cc5ccc(cc1[nH]2)n5)[nH]4)n3 '
               'NAME: 2-formylporphyrin REASON: MISSED No porphyrin ring '
               'system detected\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C(NC1=C(C(=O)CCNC(=O)CCC(NC(=O)C)C(=O)OCC(O)CO)C=CC=C1)C(C=C)(C)C',
                                     'name': 'Citrinamide B',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)CC(O)=O)C(O)=O)C',
                                     'name': 'Trp-Asp-Met',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'C(CC)(O)(C)C',
                                     'name': '2-methylbutan-2-ol',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'O1[C@@H]([C@H](O)[C@@H](O)[C@@H](NC(=O)C)C1O)CO',
                                     'name': 'N-Acetyl-D-Gulosamine',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)OS([O-])(=O)=O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)O)[H])C)[H])[H]',
                                     'name': '(24S)-hydroxycholesterol '
                                             '3-sulfate(1-)',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(C)O',
                                     'name': '17,20-dihydroxypregn-4-en-3-one',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'O=C1N2[C@H](C(=O)NCC(=O)N[C@H](C(=O)NCC(N[C@H](C(NCC(N[C@H]1CC=3C4=C(C=CC=C4)NC3)=O)=O)C(O)C5=CC(=C(O)C=C5)CC=C(C)C)=O)C(C)C)CCC2',
                                     'name': 'WIN 66306',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'CCCCCCCCCCc1c(C)c(O)c(OC)c(OC)c1O',
                                     'name': '6-decylubiquinol',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'COC1=CC(=C(C=C1Br)OC)S(=O)(=O)N2CCCCCC2',
                                     'name': '1-(4-bromo-2,5-dimethoxyphenyl)sulfonylazepane',
                                     'reason': 'No porphyrin ring system '
                                               'detected'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)N(C)C)O[C@@H]1CN(C)C(=O)NC3=CC=C(C=C3)F)[C@@H](C)CO',
                                     'name': '1-[[(2S,3R)-8-(dimethylamino)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-2-yl]methyl]-3-(4-fluorophenyl)-1-methylurea',
                                     'reason': 'No porphyrin ring system '
                                               'detected'}],
    'sample_false_negatives': [   {   'smiles': 'CC[C@@H]1[C@@H](C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(O)=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31',
                                      'name': '3-vinylbacteriochlorophyllide a',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'COC(=O)[C@H]1C(=O)c2c(C)c3C=C4C(C=C)=C(C)C5=[N+]4[Mg--]46n3c2C1=C1[C@@H](CCC(O)=O)[C@H](C)C(C=c2c(C)c(C=C)c(=C5)n42)=[N+]61',
                                      'name': 'divinyl chlorophyllide a',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'Cc1c(C=C)c2C=C3[N+]4=C(C=c5c(CCC(O)=O)c(C)c6=CC7=[N+]8C(=Cc1n2[Fe--]48n56)C(C=C)=C7C)[C@](O)(CCC(O)=O)[C@@]3(C)O',
                                      'name': 'heme d cis-diol',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'CC1OC(=O)C[C@@]2(C)[C@H](CCC(O)=O)C3=CC4=[N+]5C(=CC6=[N+]7C(=CC8=[N+]9C(=C(CC(O)=O)[C@@]8(C)CCC(O)=O)C12N3[Co--]579)C(CCC(O)=O)=C6CC(O)=O)[C@@H](CCC(O)=O)[C@]4(C)CC(O)=O',
                                      'name': 'cobalt(II)-factor IV',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'CC1=C(CCC(O)=O)C2=[N+]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N+]5[Mg--]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2',
                                      'name': 'magnesium protoporphyrin',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'Fc1c(F)c(F)c(c(F)c1F)C1=C2C=CC3=[N]2[Pd]24n5c1ccc5C(c1c(F)c(F)c(F)c(F)c1F)=C1C=CC(C(c5c(F)c(F)c(F)c(F)c5F)=c5ccc(=C3c3c(F)c(F)c(F)c(F)c3F)n25)=[N]41',
                                      'name': 'palladium(II)-5,10,15,20-tetrakis-(2,3,4,5,6-pentafluorphenyl)-porphyrin',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'CC1=C2[N+]3=C(C=C4N5C(=CC6=[N+]7C(Cc8c(CCC(O)=O)c(CC(O)=O)c1n8[Co--]357)=C(CCC(O)=O)C6CC(O)=O)[C@@H](CCC(O)=O)[C@]4(C)CC(O)=O)[C@@H](CCC(O)=O)[C@]2(C)CC(O)=O',
                                      'name': 'cobalt-precorrin-3',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': '[H]C(=O)c1c(C)c2\\C=C3/N=C([C@@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]3C)C3=c4c(C(=O)[C@@H]3C(=O)OC)c(C)c3=CC5=N\\C(=C/c1n2[Mg]n43)\\C(C)=C5CC',
                                      'name': 'chlorophyll d',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'C=1C=2C(=C(C=3N2[Pd]45N=6C(=CC=7N4C(C3)=C(C7CC)CC)C(C(C6C=C8N5C1C(=C8CC)CC)=O)(CC)CC)CC)CC',
                                      'name': 'palladium(II) '
                                              'octaethylporphyrinketone',
                                      'reason': 'No porphyrin ring system '
                                                'detected'},
                                  {   'smiles': 'Cc1c(C=C)c2C=C3[N+]4=C(C=c5c(CCC(O)=O)c(C)c6=CC7=[N+]8C(=Cc1n2[Fe--]48n56)C(C=C)=C7C)[C@@]1(CCC(=O)O1)[C@@]3(C)O',
                                      'name': 'trans-heme d hydroxychlorin '
                                              'gamma-spirolactone',
                                      'reason': 'No porphyrin ring system '
                                                'detected'}],
    'attempt': 2,
    'reasoning': 'The previous program attempted to identify porphyrins by '
                 'searching for a specific SMARTS pattern representing the '
                 'porphyrin core. However, this pattern was too specific and '
                 'did not account for variations in porphyrin structures, '
                 'leading to false negatives.\n'
                 '\n'
                 'Porphyrins are characterized by a macrocyclic structure '
                 'composed of four pyrrole rings connected through methine '
                 'bridges, forming a 16-membered ring containing four nitrogen '
                 'atoms. To improve the program, we can:\n'
                 '\n'
                 '- Search for 16-membered rings in the molecule.\n'
                 '- Check if any of these rings contain exactly four nitrogen '
                 'atoms.\n'
                 '- If such a ring is found, classify the molecule as a '
                 'porphyrin.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 55,
    'num_false_positives': 100,
    'num_true_negatives': 77317,
    'num_false_negatives': 68,
    'num_negatives': None,
    'precision': 0.3548387096774194,
    'recall': 0.44715447154471544,
    'f1': 0.39568345323741005,
    'accuracy': 0.9978333763218984,
    'negative_predictive_value': 0.999121276733217}