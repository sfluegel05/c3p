"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is any chlorobenzene carrying three chloro substituents at unspecified positions.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over all rings in the molecule
    for ring in atom_rings:
        # Check if the ring is aromatic and has 6 members
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Check if all atoms in the ring are carbons
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Check if the ring is not fused (atoms are only part of one ring)
                if all(ring_info.NumAtomRings(idx) == 1 for idx in ring):
                    # Count the number of chlorine atoms attached directly to the ring
                    num_cl = 0
                    for idx in ring:
                        atom = mol.GetAtomWithIdx(idx)
                        for neighbor in atom.GetNeighbors():
                            # Neighbor not in ring and is chlorine
                            if neighbor.GetAtomicNum() == 17 and neighbor.GetIdx() not in ring:
                                num_cl += 1
                    if num_cl == 3:
                        return True, "Contains standalone benzene ring with exactly three chlorine substituents"
    
    # No suitable benzene ring with exactly three chlorine substituents found
    return False, "Does not contain standalone benzene ring with exactly three chlorine substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27096',
                          'name': 'trichlorobenzene',
                          'definition': 'Any member of the class of '
                                        'chlorobenzenes carrying three chloro '
                                        'substituents at unspecified '
                                        'positions.',
                          'parents': ['CHEBI:23132'],
                          'xrefs': ['Wikipedia:Trichlorobenzene'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6349206349206349 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: OC(=O)C1=CC(Cl)=C(Cl)C(Cl)=C1 NAME: '
               '3,4,5-trichlorobenzoic acid REASON: CORRECT Contains benzene '
               'ring with exactly three chlorine substituents\n'
               ' * SMILES: Clc1ccc(-c2cc(Cl)c(Cl)c(Cl)c2)c(Cl)c1Cl NAME: '
               "2,3,3',4,4',5'-Hexachlorobiphenyl REASON: CORRECT Contains "
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: Oc1cc(Cl)c(cc1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "4-Hydroxy-2,2',4',5,5'-pentachlorobiphenyl REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1ccc(S\\N=N\\c2cc(Cl)c(Cl)cc2Cl)cc1 NAME: '
               'Chlorfensulphide REASON: CORRECT Contains benzene ring with '
               'exactly three chlorine substituents\n'
               ' * SMILES: C1=CC(=C(C(=C1Cl)C(O)=O)Cl)Cl NAME: 2,3,6-TBA '
               'REASON: CORRECT Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: Oc1c(Cl)cc(c(Cl)c1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "2,2',3,4',5,5'-Hexachloro-4-biphenylol REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl NAME: '
               "2,2',3,3',5,5'-hexachlorobiphenyl REASON: CORRECT Contains "
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: Oc1cc(Cl)c(Cl)cc1Cl NAME: 2,4,5-trichlorophenol '
               'REASON: CORRECT Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: COc1c(Cl)cc(Cl)cc1Cl NAME: 2,4,6-trichloroanisole '
               'REASON: CORRECT Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: Oc1c(Cl)cc(c(Cl)c1Cl)-c1ccc(Cl)c(Cl)c1Cl NAME: '
               "2,2',3,3',4',5-Hexachloro-4-biphenylol REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1ccc(Cl)c(-c2ccccc2)c1Cl NAME: '
               '2,3,6-trichlorobiphenyl REASON: CORRECT Contains benzene ring '
               'with exactly three chlorine substituents\n'
               ' * SMILES: Clc1ccc(cc1Cl)-c1cc(Cl)c(Cl)c(Cl)c1 NAME: '
               "3,3',4,4',5-pentachlorobiphenyl REASON: CORRECT Contains "
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: Clc1cc(cc(Cl)c1Cl)-c1cc(Cl)c(Cl)c(Cl)c1 NAME: '
               "3,4,5,3',4',5'-Hexachlorobiphenyl REASON: CORRECT Contains "
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: Oc1c(Cl)cc(Cl)c(Cl)c1-c1cc(Cl)cc(Cl)c1Cl NAME: '
               "2,2',3,3',5,5'-Hexachloro-6-biphenylol REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: OC(=O)C1=C(Cl)C(Cl)=CC(Cl)=C1 NAME: '
               '2,3,5-trichlorobenzoic acid REASON: CORRECT Contains benzene '
               'ring with exactly three chlorine substituents\n'
               ' * SMILES: Clc1ccc(Cl)c(Cl)c1 NAME: 1,2,4-trichlorobenzene '
               'REASON: CORRECT Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: Clc1ccccc1-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "2,2',4,5-tetrachlorobiphenyl REASON: CORRECT Contains benzene "
               'ring with exactly three chlorine substituents\n'
               ' * SMILES: Clc1cccc(Cl)c1Cl NAME: 1,2,3-trichlorobenzene '
               'REASON: CORRECT Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: Oc1c(Cl)c(Cl)c(c(Cl)c1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "4-Hydroxy-2,2',3,4',5,5',6-heptachlorobiphenyl REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Oc1ccc(c(Cl)c1)-c1c(Cl)cc(Cl)cc1Cl NAME: '
               "4-Hydroxy-2,2',4',6'-tetrachlorobiphenyl REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: OC(=O)COc1cc(Cl)c(Cl)cc1Cl NAME: '
               '(2,4,5-trichlorophenoxy)acetic acid REASON: CORRECT Contains '
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: Oc1cc(Cl)c(Cl)c(Cl)c1-n1cc(Cl)c(Cl)c1Cl NAME: '
               'neopyrrolomycin B REASON: CORRECT Contains benzene ring with '
               'exactly three chlorine substituents\n'
               ' * SMILES: Oc1c(Cl)c(Cl)cc(c1Cl)-c1cc(Cl)c(Cl)c(Cl)c1Cl NAME: '
               "2,2',3',4,4',5,5'-Heptachloro-3-biphenylol REASON: CORRECT "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1ccc(cc1)-c1c(Cl)cc(Cl)cc1Cl NAME: '
               "2,4,4',6-tetrachlorobiphenyl REASON: CORRECT Contains benzene "
               'ring with exactly three chlorine substituents\n'
               ' * SMILES: Clc1ccc(cc1Cl)-c1ccc(Cl)c(Cl)c1Cl NAME: '
               "2,3,4,3',4'-Pentachlorobiphenyl REASON: CORRECT Contains "
               'benzene ring with exactly three chlorine substituents\n'
               'False positives: SMILES: '
               'Oc1c(Cl)cc2c(oc3cc(Cl)c(Cl)c(Cl)c23)c1Cl NAME: '
               '7-hydroxy-1,2,3,6,8-pentachlorodibenzofuran REASON: WRONGLY '
               'CLASSIFIED Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: C=1(C(=CC(Cl)=CC1Cl)Cl)OC=2C=CC(=CC2)[N+]([O-])=O '
               'NAME: chlornitrofen REASON: WRONGLY CLASSIFIED Contains '
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: COP(=S)(OC)Oc1cc(Cl)c(Cl)cc1Cl NAME: Fenchlorphos '
               'REASON: WRONGLY CLASSIFIED Contains benzene ring with exactly '
               'three chlorine substituents\n'
               ' * SMILES: CC(Cl)(Cl)C(=O)OCCOc1cc(Cl)c(Cl)cc1Cl NAME: Erbon '
               'REASON: WRONGLY CLASSIFIED Contains benzene ring with exactly '
               'three chlorine substituents\n'
               ' * SMILES: Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2cc1Cl NAME: '
               '1,2,3,7,8-Pentachlorodibenzodioxin REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: C1=CN=CC=C1NC(=O)NC2=C(C=C(C=C2Cl)Cl)Cl NAME: '
               '1-pyridin-4-yl-3-(2,4,6-trichlorophenyl)urea REASON: WRONGLY '
               'CLASSIFIED Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: ClC1=C(OCCNCCC)C(Cl)=CC(=C1)Cl NAME: Prochloraz '
               'Metabolite BTS 40348 REASON: WRONGLY CLASSIFIED Contains '
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: '
               'CN1c2cc3Oc4cc5=[N+](C)C(C)(C)C=C(CS([O-])(=O)=O)c5cc4=C(c3cc2C(CS([O-])(=O)=O)=CC1(C)C)c1c(Cl)c(SCC(=O)NCCCCCC(=O)ON2C(=O)CCC2=O)c(Cl)c(Cl)c1C([O-])=O '
               'NAME: Alexa Fluor 610-X(2-) REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: '
               'C1=CC(=CC(=C1)[N+](=O)[O-])C(=O)NNC2=C(C=C(C=C2Cl)Cl)Cl NAME: '
               "3-nitro-N'-(2,4,6-trichlorophenyl)benzohydrazide REASON: "
               'WRONGLY CLASSIFIED Contains benzene ring with exactly three '
               'chlorine substituents\n'
               ' * SMILES: '
               'CCCN(CC1CC1)C2=NC(=NC(=C2Cl)NC3=C(C=C(C=C3Cl)Cl)Cl)C NAME: '
               '5-chloro-N4-(cyclopropylmethyl)-2-methyl-N4-propyl-N6-(2,4,6-trichlorophenyl)pyrimidine-4,6-diamine '
               'REASON: WRONGLY CLASSIFIED Contains benzene ring with exactly '
               'three chlorine substituents\n'
               ' * SMILES: ClC1=C(Cl)C(OC)=C(O)C(=C1OC)Cl NAME: '
               '2,4,5-Trichlor-3,6-dimethoxy-phenol REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: C1=CC(=C(C(=C1C2=NNC(=C2)C(=O)NCCO)Cl)Cl)Cl NAME: '
               'N-(2-hydroxyethyl)-3-(2,3,4-trichlorophenyl)-1H-pyrazole-5-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains benzene ring with exactly '
               'three chlorine substituents\n'
               ' * SMILES: Clc1cc2oc3cc(Cl)c(Cl)c(Cl)c3c2cc1Cl NAME: '
               '1,2,3,7,8-Pentachlorodibenzofuran REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: [O-]C(=O)C1=C(Cl)C(Cl)=CC(Cl)=C1 NAME: '
               '2,3,5-trichlorobenzoate REASON: WRONGLY CLASSIFIED Contains '
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: C=1(C(=C(C=C(C1Cl)O)Cl)O)Cl NAME: '
               '2,3,5-trichlorobenzene-1,4-diol REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: ClC=1C(C=2C(Cl)=CC(Cl)=C(Cl)C2)=CC(Cl)=C(Cl)C1Cl '
               'NAME: PCB180 REASON: WRONGLY CLASSIFIED Contains benzene ring '
               'with exactly three chlorine substituents\n'
               ' * SMILES: '
               'ClC1=C(Cl)C(OC)=C(Cl)C(=C1OC2=C(Cl)C(OC)=C(Cl)C(=C2OC3=C(Cl)C(OC)=C(Cl)C(=C3O)Cl)Cl)Cl '
               'NAME: Chlorophellin C REASON: WRONGLY CLASSIFIED Contains '
               'benzene ring with exactly three chlorine substituents\n'
               ' * SMILES: ClC1=C(Cl)C(O)=C(OC)C(=C1OC)Cl NAME: '
               '2,3,5-trichloro-4,6-dimethoxy-phenol REASON: WRONGLY '
               'CLASSIFIED Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1cc2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2c(Cl)c1Cl NAME: '
               '1,2,3,4,6,7,8-Heptachlorodibenzodioxin REASON: WRONGLY '
               'CLASSIFIED Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: ClC1=C(OCCN(C(=O)N)CCC)C(Cl)=CC(=C1)Cl NAME: BTS '
               '44595 REASON: WRONGLY CLASSIFIED Contains benzene ring with '
               'exactly three chlorine substituents\n'
               ' * SMILES: Clc1cc(Cl)c(cc1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "2,2',4,4',5,5'-hexachlorobiphenyl REASON: WRONGLY CLASSIFIED "
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: Clc1cc2c3cc(Cl)c(Cl)c(Cl)c3oc2c(Cl)c1Cl NAME: '
               '2,3,4,6,7,8-Hexachlorodibenzofuran REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: ClC1=C(OC)C(Cl)=CC(=C1O)Cl NAME: '
               '2,4,6-trichloro-3-methoxyphenol REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               ' * SMILES: CC(Oc1cc(Cl)c(Cl)cc1Cl)C(O)=O NAME: Fenoprop '
               'REASON: WRONGLY CLASSIFIED Contains benzene ring with exactly '
               'three chlorine substituents\n'
               ' * SMILES: [O-]C(=O)COc1cc(Cl)c(Cl)cc1Cl NAME: '
               '(2,4,5-trichlorophenoxy)acetate REASON: WRONGLY CLASSIFIED '
               'Contains benzene ring with exactly three chlorine '
               'substituents\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC(CCCCCCCCCCCC)C(O)CCCC(O)=O',
                                     'name': '5,6-dihydroxy stearic acid',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'O=C(O)C1=CC=C(CCC)C=C1',
                                     'name': '4-propylbenzoic acid',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'NC(=O)NCNC(N)=O',
                                     'name': 'methylenediurea',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'O=C1N2[C@H](C(=O)N3[C@@]1(O)[C@@H](OCC=C(C)C)C=4C5=C(C=C(OC)C=C5)N6C4[C@@H]3[C@@H](O)C(OO[C@@H]6C=C(C)C)(C)C)CCC2',
                                     'name': '26alpha-hydroxyfumitremorgin A',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'Cc1cccc2c(cc(O)cc12)C(O)=O',
                                     'name': '3-hydroxy-5-methyl-1-naphthoic '
                                             'acid',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'O=C(N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(C)C)C(O)=O)[C@@H](N)[C@H](CC)C',
                                     'name': 'Ile-Phe-Leu',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(18:1(9Z)/19:1(9Z))',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'P(OCC(OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'PE-NMe2(22:1(13Z)/14:0)',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12',
                                     'name': 'arachidonoyl-CoA(4-)',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H](O)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@@H]7O[C@@H]([C@@H](O[C@@H]8O[C@@H]([C@H](O)[C@H](O)[C@H]8O)CO)[C@H](O)[C@H]7NC(=O)C)CO)CO',
                                     'name': 'N-[(2S,3R,4R,5S,6R)-2-[(2R,3R,4R,5S,6S)-4-[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-Acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-5,6-dihydroxyoxan-3-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Does not contain standalone '
                                               'benzene ring with exactly '
                                               'three chlorine substituents'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'After analyzing the previous code and the outcomes, it '
                 'appears that the function `is_trichlorobenzene` is '
                 'incorrectly classifying molecules that contain benzene rings '
                 'with exactly three chlorine substituents, even when these '
                 'rings are part of fused ring systems or larger polycyclic '
                 'molecules. The main issue is that the code does not '
                 'differentiate between standalone benzene rings and benzene '
                 'rings that are fused to other rings or are part of '
                 'polycyclic systems.\n'
                 '\n'
                 'To improve the classification, we need to ensure that:\n'
                 '\n'
                 '- The benzene ring with three chlorine substituents is not '
                 'fused to any other rings (i.e., it is a standalone '
                 'monocyclic benzene ring).\n'
                 '- The atoms of the benzene ring are only part of one ring '
                 '(they are not shared with other rings).\n'
                 '- The chlorine substituents are attached directly to the '
                 'benzene ring and are not part of other functional groups or '
                 'fused systems.\n'
                 '\n'
                 'By adding a check to ensure that the benzene ring is not '
                 'fused (all its atoms are only in one ring), we can reduce '
                 'the number of false positives where the chlorinated benzene '
                 'ring is part of a larger fused system.\n'
                 '\n'
                 'Here is the updated code implementing these improvements:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 40,
    'num_false_positives': 36,
    'num_true_negatives': 142224,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.5263157894736842,
    'recall': 1.0,
    'f1': 0.6896551724137931,
    'accuracy': 0.999747013352073,
    'negative_predictive_value': 1.0}