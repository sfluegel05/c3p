"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone where the middle carbon (position 2)
    is esterified with a fatty acid, and the other two carbons (positions 1 and 3)
    have hydroxyl groups and no other substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern: three connected carbons
    glycerol_pattern = Chem.MolFromSmarts("C-C-C")
    matches = mol.GetSubstructMatches(glycerol_pattern)

    if not matches:
        return False, "No glycerol backbone found"

    # Iterate over glycerol backbone matches
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check that c1 and c3 have hydroxyl groups and no other substituents
        def is_terminal_carbon_ok(carbon_atom, connected_idx):
            oh_found = False
            for neighbor in carbon_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == connected_idx:
                    continue  # Skip the bond to the central carbon
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    # Oxygen with only one bond (hydroxyl group)
                    oh_found = True
                else:
                    return False  # Other substituent found
            return oh_found

        c1_ok = is_terminal_carbon_ok(c1, c2_idx)
        c3_ok = is_terminal_carbon_ok(c3, c2_idx)

        if not (c1_ok and c3_ok):
            continue  # Positions 1 and 3 must have only hydroxyl groups

        # Check that c2 is connected to an ester group and no other substituents
        ester_found = False
        for neighbor in c2.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in [c1_idx, c3_idx]:
                continue  # Skip bonds to the glycerol backbone carbons
            if neighbor.GetAtomicNum() == 8:
                # Check if this oxygen is part of an ester linkage
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetIdx() == c2_idx:
                        continue  # Skip the bond back to c2
                    if neighbor2.GetAtomicNum() == 6:
                        # Check for carbonyl (C=O)
                        for bond in neighbor2.GetBonds():
                            other_atom = bond.GetOtherAtom(neighbor2)
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and other_atom.GetAtomicNum() == 8:
                                ester_found = True
                                break
                        if ester_found:
                            break
                if ester_found:
                    continue
                else:
                    return False, "Position 2 oxygen is not part of ester linkage"
            else:
                return False, "Position 2 has unexpected substituents"

        if not ester_found:
            continue  # Ester linkage at position 2 is required

        # Check that there is only one ester group in the molecule
        ester_pattern = Chem.MolFromSmarts("C(=O)O")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) != 1:
            return False, f"Found {len(ester_matches)} ester groups, expected 1"

        # Ensure there are exactly two hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) != 2:
            return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, expected 2"

        # Passed all checks
        return True, "Molecule is a 2-monoglyceride with acyl group at position 2"

    return False, "Molecule does not match the 2-monoglyceride structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17389',
                          'name': '2-monoglyceride',
                          'definition': 'A monoglyceride in which the acyl '
                                        'substituent is located at position 2.',
                          'parents': ['CHEBI:17408', 'CHEBI:76575'],
                          'xrefs': ['KEGG:C02112'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.544 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CCCCCCCCCCCC(=O)OC(CO)CO NAME: '
               '2-lauroylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/CCCCCCCC(=O)OC(CO)CO NAME: '
               '2-alpha-linolenoylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0) REASON: CORRECT '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C\\C=C/C=C/[C@H](CCCCC)OO)OC(CO)CO NAME: '
               '15(S)-HPETE 2-glyceryl ester REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/22:1(13Z)/0:0) REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCC=CCC=CCCCCC)C(CO)CO NAME: '
               '2-Linoleoyl glycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: CCCCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO NAME: '
               '2-oleoylglycerol REASON: CORRECT Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               ' * SMILES: CCCCCCCCCCCCCC(=O)OC(CO)CO NAME: '
               '2-tetradecanoylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC NAME: '
               'MG(0:0/22:5(4Z,7Z,10Z,13Z,16Z)/0:0) REASON: CORRECT Molecule '
               'is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'C(\\[C@H]1[C@@H](CC([C@@H]1C/C=C\\CCCC(OC(CO)CO)=O)=O)O)=C/[C@H](CCCCC)O '
               'NAME: prostaglandin E2 2-glyceryl ester REASON: CORRECT '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/15:0/0:0) REASON: CORRECT Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)OC(CO)CO NAME: '
               '2-[(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl]-sn-glycerol REASON: '
               'CORRECT Molecule is a 2-monoglyceride with acyl group at '
               'position 2\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(=O)OC(CO)CO NAME: '
               '2-palmitoylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO NAME: '
               'MG(0:0/22:2(13Z,16Z)/0:0) REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'CCCCC[C@H](OO)\\C=C\\[C@H]1[C@H]2C[C@H](OO2)[C@@H]1C\\C=C/CCCC(=O)OC(CO)CO '
               'NAME: prostaglandin G2 2-glyceryl ester REASON: CORRECT '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCC/C=C\\CCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/24:1(15Z)/0:0) REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O(C(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)C(CO)CO '
               'NAME: MG(0:0/24:6(6Z,9Z,12Z,15Z,18Z,21Z)/0:0) REASON: CORRECT '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)OC(CO)CO NAME: '
               '2-linoleoylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: OCC(OC(=O)CCCCCCC/C=C/CCCCCCCC)CO NAME: '
               '2-[(9E)-9-octadecenoyl]glycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'OCC(CO)OC(CCC/C=C/1\\O[C@]2(C[C@H]([C@@H]([C@]2(C1)[H])/C=C/[C@H](CCCCC)O)O)[H])=O '
               'NAME: prostaglandin I2 2-glyceryl ester REASON: CORRECT '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO '
               'NAME: 12-HPETE 2-glyceryl ester REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: CCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO NAME: '
               '2-[(9Z)-hexadecenoyl]glycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\CCCCCC)C(CO)CO NAME: '
               'MG(0:0/18:1(11Z)/0:0) REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               '[H]C(CO)(CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC NAME: '
               '2-arachidonoylglycerol REASON: CORRECT Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/22:0/0:0) REASON: CORRECT Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               'False positives: SMILES: '
               'O[C@@H]1[C@H]([C@@H](C(=O)C1)C/C=C\\CCCC(OC(CO)CO)=O)/C=C/[C@@H](O)CCCCC '
               'NAME: Prostaglandin PGE2 glyceryl ester REASON: WRONGLY '
               'CLASSIFIED Molecule is a 2-monoglyceride with acyl group at '
               'position 2\n'
               ' * SMILES: '
               'COc1cc(cc(OC)c1O)\\C=C\\C(=O)O[C@H]([C@@H](O)[C@H](O)[C@H](O)C([O-])=O)C([O-])=O '
               'NAME: 2-O-sinapoyl-D-glucarate(2-) REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O[C@@H]([C@H](O)[C@@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)[C@H](O)C(O)=O '
               'NAME: 2-O-caffeoylglucaric acid REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'C(=O)([C@H](OC(C=1C=C(C(=C(C1)O)O)O)=O)[C@@H](O)[C@H](O)[C@H](O)CO)[H] '
               'NAME: 2-galloyl-D-glucose REASON: WRONGLY CLASSIFIED Molecule '
               'is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)C(O)=O NAME: '
               '(2S,3R)-trans-caftaric acid REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: CCCCCC=CCC=CCC=CCC=CCCCC(=O)OC(CO)CO NAME: '
               'eicosa-5,8,11,14-tetraenoic acid 1,3-dihydroxypropan-2-yl '
               'ester REASON: WRONGLY CLASSIFIED Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:4(8Z,11Z,14Z,17Z)/0:0) REASON: WRONGLY '
               'CLASSIFIED Molecule is a 2-monoglyceride with acyl group at '
               'position 2\n'
               ' * SMILES: '
               'O([C@]([C@@H](O)CO)(C)C(O)=O)C(=O)/C=C/C1=CC(O)=C(O)C=C1 NAME: '
               '(2s,3s)-2-{[(2e)-3-(3,4-Dihydroxyphenyl)-2-propenoyl]oxy}-3,4-dihydroxy-2-methylbutanoic '
               'acid REASON: WRONGLY CLASSIFIED Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C '
               'NAME: Roselipin 3E REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: COc1ccc(CC(O)C(OC(C)=O)C(O)CO)cc1 NAME: karalicin '
               'REASON: WRONGLY CLASSIFIED Molecule is a 2-monoglyceride with '
               'acyl group at position 2\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H]([C@H]([C@H](O)C1)/C=C/[C@@H](O)CCCCC)CC(=O)CCCCC(OC(CO)CO)=O '
               'NAME: 2-glyceryl-6-keto-PGF1alpha REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O=C(OC(CO)CO)C1=C(O)C=C(O)C=C1C NAME: '
               '2-(2,4-dihydroxy-6-methylbenzoyl)-glycerol REASON: WRONGLY '
               'CLASSIFIED Molecule is a 2-monoglyceride with acyl group at '
               'position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\CCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/20:1(11Z)/0:0) REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: OC(C(OC(=O)C=Cc1ccc(O)cc1)C(O)=O)C(O)=O NAME: '
               'p-coutaric acid REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O[C@](CC1=CC(O)=C(O)C=C1)([C@H](OC(=O)/C=C/C2=CC(O)=C(O)C=C2)C(O)=O)C(O)=O '
               'NAME: Fukinolic acid REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O(C(C(O)(C)C)C(O)C1=C(OC)C=C2OC(=O)C=CC2=C1)C(=O)C=C(C)C NAME: '
               '[1,3-Dihydroxy-1-(7-methoxy-2-oxochromen-6-yl)-3-methylbutan-2-yl] '
               '3-methylbut-2-enoate REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'OC(CCCCCC/C=C/[C@H](O)[C@@H](OC(=O)C)[C@H](O)[C@H](N)C(O)=O)CCCCCC '
               'NAME: '
               '(E,2S,3R,4R,5S)-4-Acetyloxy-2-amino-3,5,14-trihydroxyicos-6-enoic '
               'acid REASON: WRONGLY CLASSIFIED Molecule is a 2-monoglyceride '
               'with acyl group at position 2\n'
               ' * SMILES: COc1cc(cc(OC)c1O)\\C=C\\C(=O)OC(C(O)=O)C(O)=O NAME: '
               'sinapoyltartronic acid REASON: WRONGLY CLASSIFIED Molecule is '
               'a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O[C@H]([C@H](OC(=O)\\C=C/c1ccc(O)c(O)c1)C(O)=O)C(O)=O NAME: '
               '(2S,3R)-cis-caftaric acid REASON: WRONGLY CLASSIFIED Molecule '
               'is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'C1=CC(=C(C=C1/C=C/C(OC(C(C(C)C)(C(=O)O)O)C(O)=O)=O)O)O NAME: '
               '3-O-caffeoyl-2-(2-propyl)-2-hydroxybutanedioic acid REASON: '
               'WRONGLY CLASSIFIED Molecule is a 2-monoglyceride with acyl '
               'group at position 2\n'
               ' * SMILES: '
               'O[C@@H](CCCCC)/C=C/C=C\\C/C=C\\C/C=C\\CCCC(OC(CO)CO)=O NAME: '
               '15-HETE-G REASON: WRONGLY CLASSIFIED Molecule is a '
               '2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)cc1)C(O)=O)C(O)=O '
               'NAME: (2R,3S)-trans-coutaric acid REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO NAME: '
               'MG(0:0/20:2(11Z,14Z)/0:0) REASON: WRONGLY CLASSIFIED Molecule '
               'is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:3(11Z,14Z,17Z)/0:0) REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               ' * SMILES: '
               'O[C@@H]([C@H](O)[C@@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C([O-])=O)[C@H](O)C([O-])=O '
               'NAME: 2-O-caffeoylglucarate(2-) REASON: WRONGLY CLASSIFIED '
               'Molecule is a 2-monoglyceride with acyl group at position 2\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Molecule does not match the '
                                               '2-monoglyceride structure'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(CCC/C=C\\C/C=C\\C\\C=C/C=C/[C@H](CCCCC)OO)OC(CO)CO',
                                      'name': '15(S)-HPETE 2-glyceryl ester',
                                      'reason': 'Found 3 hydroxyl groups, '
                                                'expected 2'},
                                  {   'smiles': 'C(\\[C@H]1[C@@H](CC([C@@H]1C/C=C\\CCCC(OC(CO)CO)=O)=O)O)=C/[C@H](CCCCC)O',
                                      'name': 'prostaglandin E2 2-glyceryl '
                                              'ester',
                                      'reason': 'Found 4 hydroxyl groups, '
                                                'expected 2'},
                                  {   'smiles': 'CCCCC[C@H](OO)\\C=C\\[C@H]1[C@H]2C[C@H](OO2)[C@@H]1C\\C=C/CCCC(=O)OC(CO)CO',
                                      'name': 'prostaglandin G2 2-glyceryl '
                                              'ester',
                                      'reason': 'Found 3 hydroxyl groups, '
                                                'expected 2'},
                                  {   'smiles': 'OCC(CO)OC(CCC/C=C/1\\O[C@]2(C[C@H]([C@@H]([C@]2(C1)[H])/C=C/[C@H](CCCCC)O)O)[H])=O',
                                      'name': 'prostaglandin I2 2-glyceryl '
                                              'ester',
                                      'reason': 'Found 4 hydroxyl groups, '
                                                'expected 2'},
                                  {   'smiles': 'O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO',
                                      'name': '12-HPETE 2-glyceryl ester',
                                      'reason': 'Found 3 hydroxyl groups, '
                                                'expected 2'},
                                  {   'smiles': 'C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO',
                                      'name': 'prostaglandin D2 2-glyceryl '
                                              'ester',
                                      'reason': 'Found 4 hydroxyl groups, '
                                                'expected 2'}],
    'attempt': 1,
    'reasoning': 'Upon analyzing the previous program and the outcomes, it '
                 'appears that the code incorrectly classified certain '
                 'molecules as 2-monoglycerides. The main issues seem to be:\n'
                 '\n'
                 '1. **Overly Broad SMARTS Patterns**: The SMARTS patterns '
                 'used were not specific enough, leading to false positives. '
                 'Molecules with similar substructures or additional '
                 'functional groups were incorrectly matched.\n'
                 '\n'
                 '2. **Insufficient Checks for Additional Functional Groups**: '
                 'The code did not adequately check for other ester or acyl '
                 'groups elsewhere in the molecule, resulting in '
                 'misclassification of molecules that are not true '
                 '2-monoglycerides.\n'
                 '\n'
                 '3. **Inaccurate Identification of Positions**: The code may '
                 'not have correctly identified the ester linkage at position '
                 '2 and hydroxyl groups at positions 1 and 3.\n'
                 '\n'
                 'To improve the classification, we can:\n'
                 '\n'
                 '- **Use Specific SMARTS Patterns**: Define more precise '
                 'SMARTS patterns that accurately represent the structure of '
                 '2-monoglycerides, ensuring that the glycerol backbone and '
                 'functional groups are correctly identified.\n'
                 '\n'
                 '- **Check for Additional Functional Groups**: Include checks '
                 'to ensure there are no extra ester or acyl groups in the '
                 'molecule that would disqualify it from being a '
                 '2-monoglyceride.\n'
                 '\n'
                 '- **Verify Connectivity and Substituents**: Ensure that '
                 'positions 1 and 3 of the glycerol backbone have hydroxyl '
                 'groups and no other substituents, and that position 2 has an '
                 'ester linkage to a fatty acid chain.\n'
                 '\n'
                 'Here is the improved code incorporating these changes:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 28,
    'num_false_positives': 6,
    'num_true_negatives': 142260,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.8235294117647058,
    'recall': 0.8235294117647058,
    'f1': 0.8235294117647058,
    'accuracy': 0.9999156711173577,
    'negative_predictive_value': 0.9999578254818439}