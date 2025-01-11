"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all rings in the molecule
    ssr = Chem.GetSymmSSSR(mol)
    macrocycle_found = False

    # Iterate over each ring
    for ring in ssr:
        ring_size = len(ring)
        # Check for rings with 12 or more atoms
        if ring_size >= 12:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_atom_idxs = set(ring)
            # Search for lactone functionality within the ring
            for atom in ring_atoms:
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check if carbon is a carbonyl carbon within the ring
                    is_carbonyl = False
                    for nbr in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                            is_carbonyl = True  # Found C=O
                            break
                    if is_carbonyl:
                        # Check for single bond to oxygen within the ring (ester linkage)
                        for nbr in atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nbr.GetAtomicNum() == 8:
                                # Ensure the oxygen is in the ring
                                if nbr.GetIdx() in ring_atom_idxs:
                                    macrocycle_found = True
                                    return True, "Contains a macrocyclic lactone ring of 12 or more members"
    if not macrocycle_found:
        return False, "No macrocyclic lactone ring of 12 or more members found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25106',
                          'name': 'macrolide',
                          'definition': 'A macrocyclic lactone with a ring of '
                                        'twelve or more members derived from a '
                                        'polyketide.',
                          'parents': ['CHEBI:26188', 'CHEBI:63944'],
                          'xrefs': ['Wikipedia:Macrolide'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.3902439024390244 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCN1C[C@H](C)[C@@H](O)[C@](C)(O)[C@@H](CC)OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@H]1C '
               'NAME: gamithromycin REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               '[C@@]12(O[C@@H]([C@@H](C)C=C1)C3CCCCC3)C[C@@H]4C[C@](O2)(CC=C([C@@H](O[C@H]5C[C@@H]([C@@H](O[C@H]6CC([C@@H](O)[C@@H](O6)C)OC)[C@@H](O5)C)OC)[C@H](C=CC=C7[C@]8([C@H](C(O4)=O)C=C(C)[C@H]([C@]8(OC7)[H])O)O)C)C)[H] '
               'NAME: Doramectin REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'O1C2(OC3CC(OC(=O)C4C(O)(CC(OC)C(=C4)C)C(CO)=CC=CC(CC(=CC3)C)C)C2)CCC(C1CC)C '
               'NAME: Milbemycin beta2 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'S1C(=NC(=C1)/C=C(/[C@H]2OC(=O)C[C@H](O)C(C(=O)[C@@H]([C@@H](O)[C@@H](C)CC(CCCC(=CC2)C)C)C)(C)C)\\C)C '
               'NAME: Epothilone I3 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: S(C1C(=O)C(O)CCCCCC(OC(C1)=O)C)CC(O)C(=O)OC NAME: '
               'Pandangolide 3 REASON: CORRECT Contains a macrocyclic lactone '
               'ring of 12 or more members\n'
               ' * SMILES: '
               '[H][C@@]1(OC(=O)\\C=C\\CC\\C=C\\[C@H](OC)[C@@H](O)[C@H](C)\\C=C1\\C)[C@H](C)C(=O)CCCC1CC(=O)NC(=O)C1 '
               'NAME: migrastatin REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'S1C(=NC(=C1)/C=C(/[C@H]2OC(=O)C[C@H](O)C(C(=O)[C@@H]([C@@H](O)[C@H](CCCC=CC2)C)C)C)\\C)C '
               'NAME: Epothilone C1 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               '[C@@H]1([C@H]([C@@]2([C@@H](OC(=O)C)[C@H]([C@]3([C@](O)([C@]42[C@H](OC(=O)C)[C@]1([H])[C@@](COC(C5=C(CC[C@@H](C(O3)=O)C)N=CC=C5)=O)(C)O4)C)[H])O)COC(=O)C)OC(=O)C)OC(C)=O '
               'NAME: 2-O-deacetyleuonine REASON: CORRECT Contains a '
               'macrocyclic lactone ring of 12 or more members\n'
               ' * SMILES: '
               'S(C1C(=O)C(O)CCCCCC(OC(C1)=O)C)C2C(=O)C(O)CCCCCC(OC(C2)=O)C '
               'NAME: Pandangolide 4 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC(C)=O)[C@@H](OC(=O)CC)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@@](C)(C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O)O[C@H]1C[C@H]([C@@H](O)[C@H](C)O1)N(C)C '
               'NAME: Megalomicin C2 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'S1C(=NC(=C1)/C=C(/[C@H]2OC(=O)C[C@H](O)C(C(=O)[C@@H]([C@@H](O)CCCCC=CC2)C)(C)C)\\C)C '
               'NAME: Epothilone C4 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]/2OC=C3[C@](\\C2=C\\C)(CC(OCCC4=CC=C(OC=5C=C(CCOC3=O)C=CC5O)C=C4)=O)[H])CO '
               'NAME: '
               '(15S,16Z,17S)-16-ethylidene-4-hydroxy-15-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,10,14,20-tetraoxatetracyclo[21.2.2.13,7.012,17]octacosa-1(25),3,5,7(28),12,23,26-heptaene-11,19-dione '
               'REASON: CORRECT Contains a macrocyclic lactone ring of 12 or '
               'more members\n'
               ' * SMILES: '
               'O=C1C=CC(=C[C@H]([C@H](OC(=O)C[C@H]([C@@H]([C@H]([C@H](C[C@H]1C)CC)O[C@@H]2O[C@H]([C@H](O)[C@@](C2)(O)C)C)C)O)CC)C)C '
               'NAME: 5-O-mycarosyltylactone REASON: CORRECT Contains a '
               'macrocyclic lactone ring of 12 or more members\n'
               ' * SMILES: '
               'CCCCCC(=O)O[C@H]1[C@H](C)O[C@H](C[C@@]1(C)O)O[C@@H]1[C@@H](C)O[C@@H](O[C@H]2[C@@H](CC=O)C[C@@H](C)[C@@H](O)\\C=C\\C=C\\C[C@@H](C)OC(=O)C[C@@H](O)[C@@H]2OC)[C@H](O)[C@H]1N(C)C '
               'NAME: Leucomycin A13 REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'CCC(O)CC1CCCC2(CC3OC(=O)\\C=C\\C(C)(O)C(O)C(C)C(O)C(O)C(O)C(C)(O)CCCCC\\C=C\\C4CC(C)(C)OC4(O)CC(O2)C3C)O1 '
               'NAME: NK154183A REASON: CORRECT Contains a macrocyclic lactone '
               'ring of 12 or more members\n'
               ' * SMILES: '
               'CC[C@@H](O)C[C@H]1CCC[C@@]2(C[C@@H]3OC(=O)\\C=C\\[C@@](C)(O)[C@@H](O)[C@H](C)[C@@H](O)[C@H](O[C@H]4CC[C@@H]([C@H](C)O4)N(C)C)[C@@H](O)[C@](C)(O)CCCCC\\C=C\\[C@@H]4CC(C)(C)O[C@@]4(O)C[C@H](O2)[C@H]3C)O1 '
               'NAME: ossamycin REASON: CORRECT Contains a macrocyclic lactone '
               'ring of 12 or more members\n'
               ' * SMILES: '
               '[H][C@@]1(OC(=O)\\C(OC)=C\\C(C)=C\\[C@@H](C)[C@@H](O)[C@@H](CC)[C@@H](O)[C@H](C)C\\C(C)=C\\C=C\\[C@@H]1OC)[C@@H](C)[C@@H](O)[C@H](C)[C@@]1(O)C[C@@H](O[C@H]2C[C@@H](O)[C@H](OC(N)=O)[C@@H](C)O2)[C@H](C)[C@H](O1)\\C=C\\C '
               'NAME: concanamycin A REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'ClC1=C2C3=CC(=C1)[C@@H](O)COC(=O)[C@H](OC)[C@]45[C@H]([C@@]6(O[C@]37[C@H](O)C[C@H]([C@@]7([C@@H]2O)OC6=C(C4=O)C)O)O)O5 '
               'NAME: Sporolide B REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'O1[C@]2(O[C@@H]([C@H](C=C2)C)C(C)C)C[C@]3(OC(=O)[C@]4([C@@]5(O)[C@](OCC5=CC=C[C@@H]([C@H](O)C(=CC[C@@]1(C3)[H])C)C)(C(=O)C(=C4)C)[H])[H])[H] '
               "NAME: 5-Oxoavermectin ''1b'' aglycone REASON: CORRECT Contains "
               'a macrocyclic lactone ring of 12 or more members\n'
               ' * SMILES: '
               'BrC=C(OC)C(=C)CC(O)CC(O)C(C(O)CC(O)C(C(O)CC(=C[C@H](O)[C@]1(O[C@H]2CC(CC(O)CC(CC(=CC(O[C@@H]1C2)=O)C)C)C)C)C)(C)C)C '
               'NAME: Oscillariolide REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'ClC1=C(O)C=C(OC)C2=C1CC(=O)C[C@@H](OC)CCC(CC[C@@H](OC2=O)C)=O '
               'NAME: Rhinoclactone D REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               '[H][C@]12C[C@@H](O[C@@H]3O[C@H](C)[C@@H](O)[C@H](N)[C@@H]3O)\\C=C\\C=C\\C=C\\C=C\\C=C\\C=C\\C=C\\[C@H](C)[C@@H](O)[C@@H](C)[C@H](C)OC(=O)C[C@H](O)C[C@H](O)CC[C@@H](O)[C@H](O)C[C@H](O)C[C@](O)(C[C@H](O)[C@H]1C(O)=O)O2 '
               'NAME: amphotericin B REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'C\\C=C1\\C[C@@H](C)[C@@](C)(O)C(=O)OC[C@H]2CCN3CC[C@@H](OC1=O)[C@@H]23 '
               'NAME: Platyphylline REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'ClC1=C(O)C=C(O)C2=C1CC(=O)C[C@@H](O)C=CCCC=CC[C@H](OC2=O)C '
               'NAME: Ilyoresorcy B REASON: CORRECT Contains a macrocyclic '
               'lactone ring of 12 or more members\n'
               ' * SMILES: '
               'CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@H](OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@H]1C '
               'NAME: 3-O-(alpha-L-oleandrosyl)erythronolide B REASON: CORRECT '
               'Contains a macrocyclic lactone ring of 12 or more members\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C[C@@H]1CCC[C@H](C2CNC(N)=N2)C(=O)OC[C@H](C)CCC[C@H](C2CNC(N)=N2)C(=O)OC1 '
               'NAME: Chaksine REASON: MISSED No macrocyclic lactone ring of '
               '12 or more members found\n'
               ' * SMILES: '
               'O=C1O[C@H](C(=O)N([C@H]([C@H](O)CC(O[C@@H]([C@@H]1C)CCCCCCC)=O)CC2=CC=CC=C2)C)C '
               'NAME: Hapalosin G REASON: MISSED No macrocyclic lactone ring '
               'of 12 or more members found\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H]([C@@H]([C@@H](C=C(C=C1C)C)C)O)C)C)[C@@H]([C@H](O)[C@H]([C@]2(O[C@H]([C@@H](C)[C@H](C2)OC(=O)/C=C/C(=O)NC3=C(O)CCC3=O)C(C)C)O)C)C '
               'NAME: Makinolide REASON: MISSED No macrocyclic lactone ring of '
               '12 or more members found\n'
               ' * SMILES: O=C1O[C@@H](CC=CCC[C@H](CCC1)O)C NAME: '
               '5S-hydroxyrecifeiolide REASON: MISSED No macrocyclic lactone '
               'ring of 12 or more members found\n'
               ' * SMILES: O1C(C(O)(C)C=CC(=O)C(C=C(C(=O)C(CC(C1=O)C)C)C)C)CC '
               'NAME: '
               '(7E,11E)-14-ethyl-13-hydroxy-3,5,7,9,13-pentamethyl-1-oxacyclotetradeca-7,11-diene-2,6,10-trione '
               'REASON: MISSED No macrocyclic lactone ring of 12 or more '
               'members found\n'
               ' * SMILES: '
               'O=C1O[C@@H]2[C@@H]([C@@H](O[C@]3(C2)O[C@H]([C@@H](C)CC3)C[C@H](O)C)CC[C@H](C=CC=CCC[C@H]([C@H]([C@H](C[C@H](C[C@H]([C@@H]([C@@H]([C@H](C=C1)C)O)C)O)OC)O)C)O)CC)CCC '
               'NAME: Neomaclafungin F REASON: MISSED No macrocyclic lactone '
               'ring of 12 or more members found\n'
               ' * SMILES: '
               'O=C1OC(CCC(C(O)C(C=C(CC(C=C(COC(CC(C1)C(=O)O)=O)CCC(=O)C(CC(C(=O)O)C)CC)C)C)C)C)CCCC/C=C/C(CCC(OC(=O)CC(C(=O)O)CC(=O)O)C)C '
               'NAME: Actinoplanic acid A REASON: MISSED No macrocyclic '
               'lactone ring of 12 or more members found\n'
               ' * SMILES: '
               'O=C1O[C@@H](C[C@@H]2O[C@H]([C@@H](C(=O)O[C@@H](C[C@@H]3O[C@H]([C@@H](C(O[C@@H](C[C@@H]4O[C@H]([C@@H]1C)CC4)C)=O)C)CC3)C)C)CC2)CC '
               'NAME: '
               '(1S,2S,5R,7R,10S,11S,14R,16R,19S,20S,23R,25R)-5-ethyl-2,11,14,20,23-pentamethyl-4,13,22,28,29,30-hexaoxatetracyclo[23.2.1.17,10.116,19]triacontane-3,12,21-trione '
               'REASON: MISSED No macrocyclic lactone ring of 12 or more '
               'members found\n'
               ' * SMILES: '
               'O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)[C@H]5OC4)CO)O)C)C)O[C@@H]([C@@H](C)CC3)C '
               'NAME: Milbemycin alpha26 REASON: MISSED No macrocyclic lactone '
               'ring of 12 or more members found\n'
               ' * SMILES: '
               'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C=C/C)[C@@H](C)[C@@H](C2)O)O)C)C '
               'NAME: Micromonospolide B REASON: MISSED No macrocyclic lactone '
               'ring of 12 or more members found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'},
                                 {   'smiles': 'O=C1C2=C(C(OC)=CC=C2OC3=C1C(OC)=CC(=C3)C(=O)OC)C(=O)OC',
                                     'name': '2,8-dimethoxy-1,6-dimethoxycarbonyl-xanthen-9-one',
                                     'reason': 'No macrocyclic lactone ring of '
                                               '12 or more members found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1O[C@@]23[C@@](C=C(C=CCCC=C4COC(=C1C2=O)C=C4)C)(C=C(C)[C@H](C3)C)C',
                                      'name': 'Spirohexenolide B',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'O=C1O[C@@H]2C[C@@H](O)[C@@H]3O[C@H](CC=CC=C[C@@H](C[C@H]4O[C@@H](C[C@@H]([C@]1(O)C2)OC(=O)CCC)C[C@@H](O[C@@H]5O[C@H]([C@@H](O)[C@H]([C@H]5OC)O)C)C4)C)[C@@H](C3)C',
                                      'name': 'Mandelalide B',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'O=C1O[C@]23[C@H](C=C(C)[C@H](C2)C)CCC[C@@H](C)C[C@H](C=C(C(=C1C3=O)O)C)C',
                                      'name': 'Okilactomycin D',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': '[H][C@]12O[C@]1([H])C(C)(C)[C@]1([H])OC(=O)CC[C@](O)([C@H](OC(=O)C(C)C)[C@H]1OC(C)=O)[C@H](OC(C)=O)[C@@]1([H])[C@@H](OC(C)=O)[C@@](C)(C[C@]1(OC(C)=O)C(=O)[C@@H]2C)OC(C)=O',
                                      'name': 'Salicinolide',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'ClC=1NC(C(=O)N[C@H]2[C@H](O)C[C@H](O[C@@H]3[C@H]4C=C[C@H]5CCCCC=C[C@@]6([C@]7(C(C(=C([C@]5([C@H]4CC[C@@H]3C)CC)O)C(=O)O7)=O)C[C@H](CC)C(=C6)C(=O)O)C)O[C@@H]2C)=C(C1)Cl',
                                      'name': 'Nai414-A',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': '[H][C@]12CC[C@@](C)(OC(=O)C1=C)[C@@H](O)CCC(C)=CCC[C@](C)(O)[C@H](O)C2',
                                      'name': 'Capillolide',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'S1C=2C(=O)NCC(OC)C=CC=CC[C@H](CC3=NC(C(N[C@H](C1)C(=O)O[C@@H](CC2)C)=O)=CO3)O',
                                      'name': 'Griseoviridin B',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': '[H]C(=O)Nc1cccc(C(=O)NC2C(C)OC(=O)C(CCCC)C(OC(=O)Cc3ccccc3)C(C)OC2=O)c1O',
                                      'name': 'antimycin A9',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'O=C1O[C@@H]2C[C@@H](O)[C@@H]3O[C@H](CC=CC=C[C@@H](C[C@H]4O[C@@H](C[C@@H]([C@]1(O)[C@H]2OC(=O)CCC)OC(=O)CCC)C[C@@H](O)C4)C)[C@@H](C3)C',
                                      'name': 'Mandelalide D',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'},
                                  {   'smiles': 'O=C1O[C@@H]2C[C@@H](O)[C@@H]3O[C@H](CC=CC=C[C@@H](C[C@H]4O[C@@H](C[C@@H]([C@]1(O)[C@H]2O)OC(=O)C)C[C@@H](O)C4)C)[C@@H](C3)C',
                                      'name': 'Mandelalide C',
                                      'reason': 'No macrocyclic lactone ring '
                                                'of 12 or more members found'}],
    'attempt': 1,
    'reasoning': 'Analyzing the previous program and the outcomes, it appears '
                 'that the function is not correctly identifying all '
                 'macrolides. Specifically, many known macrolides are being '
                 'missed (false negatives). \n'
                 '\n'
                 '**Possible reasons for the failure:**\n'
                 '\n'
                 '1. **Lactone Pattern Matching:** The SMARTS pattern used to '
                 'identify the lactone group (`"C(=O)O"`) may not be adequate '
                 'for matching lactone functionalities within macrocyclic '
                 'rings. This pattern might miss lactones where the ester '
                 'linkage is part of a complex ring system or has '
                 'substituents.\n'
                 '\n'
                 '2. **Ring Substructure Extraction:** The method used to '
                 'extract the ring substructure (`Chem.PathToSubmol`) may not '
                 'include all necessary atoms and bonds for accurate pattern '
                 'matching, especially in complex macrocycles.\n'
                 '\n'
                 '3. **Bond Detection within Rings:** The program may not '
                 'correctly identify ester bonds within the ring, especially '
                 'if the lactone oxygen or carbonyl carbon is shared with '
                 'substituents or other rings.\n'
                 '\n'
                 '**Improvements to the Program:**\n'
                 '\n'
                 'To address these issues, we can implement the following '
                 'changes:\n'
                 '\n'
                 '- **Enhanced Lactone Detection:** Instead of relying on a '
                 'simple SMARTS pattern, we can programmatically search for '
                 'ester bonds within the ring. Specifically, we can look for a '
                 'carbonyl carbon atom (`C=O`) connected via a single bond to '
                 'an oxygen atom within the ring.\n'
                 '\n'
                 '- **Accurate Ring Atom and Bond Identification:** Ensure '
                 "that we're correctly identifying all atoms and bonds that "
                 'are part of each ring, including any heteroatoms like '
                 'oxygen.\n'
                 '\n'
                 '- **Comprehensive Ring Analysis:** Analyze each ring of size '
                 '≥12 thoroughly to detect the presence of the lactone ester '
                 'linkage, considering all possible bonding scenarios within '
                 'macrocycles.\n'
                 '\n'
                 'By implementing these improvements, we should enhance the '
                 'accuracy of macrolide classification and reduce the number '
                 'of false negatives.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1202,
    'num_false_positives': 100,
    'num_true_negatives': 6260,
    'num_false_negatives': 25,
    'num_negatives': None,
    'precision': 0.9231950844854071,
    'recall': 0.9796251018744906,
    'f1': 0.9505733491498616,
    'accuracy': 0.9835244497166206,
    'negative_predictive_value': 0.9960222752585521}