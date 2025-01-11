"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are large polyphenolic compounds composed of phenolic acids
    or catechins, often forming complex structures through ester or
    glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for multiple aromatic rings (common in tannins: catechin or ellagic units)
    aromatic_ring_pattern = Chem.MolFromSmarts('c1ccccc1')
    aromatic_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if len(aromatic_matches) < 2:
        return False, f"Found only {len(aromatic_matches)} aromatic rings, need at least 2"

    # Check for a variety of linking bond types to assert tannin-like complexity
    linkage_patterns = [
        Chem.MolFromSmarts('C(=O)O'),                              # Ester bond
        Chem.MolFromSmarts('c1cc(O)c(O)cc1'),                      # phenolic rings with hydroxyl groups
        Chem.MolFromSmarts('O[C@H]1[C@H](O)C[C@H](O)[C@@H](O)[C@@H]1O'),  # Glycosidic linkage
        Chem.MolFromSmarts('c-c(=O)-o'),                          # Approximation of ester in a ring
    ]

    # Flag for recognized linkages
    linkage_found = any(mol.HasSubstructMatch(pattern) for pattern in linkage_patterns)
    
    if not linkage_found:
        return False, "Lacking typical tannin linkages or catechin hinges"

    return True, "Molecule matches tannin profile with multiple aromatic rings and appropriate linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26848',
                          'name': 'tannin',
                          'definition': 'Any of a group of astringent '
                                        'polyphenolic vegetable principles or '
                                        'compounds, chiefly complex glucosides '
                                        'of catechol and pyrogallol.',
                          'parents': ['CHEBI:26195'],
                          'xrefs': ['Wikipedia:Tannin'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.4584527220630372 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C2C(OC(=O)C=3C(C=4C(C(OC2)=O)=CC(O)=C(O)C4O)=C(O)C(O)=C(O)C3)[C@H](O)[C@@H](O)[C@@H]1OC5=CC(O)=C(C(O)=C5)C(=O)C6=CC=CC=C6 '
               'NAME: Guavin B REASON: CORRECT Molecule matches tannin profile '
               'with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'S(OC=1C(O)=CC(C(OC2C(OC(=O)C3=CC(O)=C(O)C(O)=C3)C(OC(O)C2O)CO)=O)=CC1O)(O)(=O)=O '
               'NAME: '
               '[4-({[2,3-dihydroxy-6-(hydroxymethyl)-5-(3,4,5-trihydroxybenzoyloxy)oxan-4-yl]oxy}carbonyl)-2,6-dihydroxyphenyl]oxidanesulfonic '
               'acid REASON: CORRECT Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'Oc1cc(cc(O)c1O)C(=O)Oc1cc(cc(O)c1O)C(=O)OC[C@H]1OC(OC(=O)c2cc(O)c(O)c(OC(=O)c3cc(O)c(O)c(O)c3)c2)[C@H](OC(=O)c2cc(O)c(O)c(OC(=O)c3cc(O)c(O)c(O)c3)c2)[C@@H](OC(=O)c2cc(O)c(O)c(OC(=O)c3cc(O)c(O)c(O)c3)c2)[C@@H]1OC(=O)c1cc(O)c(O)c(OC(=O)c2cc(O)c(O)c(O)c2)c1 '
               'NAME: tannic acid REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O1C2C(O)(C(CC=3C(C=4C2=CC(OC)=C(OC)C4OC)=C(OCC(C(O)(C1=O)C)C)C=5OCOC5C3)C)C '
               'NAME: Gomisin D REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O[C@@H]1[C@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(Oc4c(O)c(O)c(O)cc4C(O)=O)cc3C(=O)O[C@@H]1[C@@H](O)[C@H](OC(=O)c1cc(O)c(O)c(O)c1)O2 '
               'NAME: Mallotinic acid REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](OC(=O)C2=CC(OC)=C(O)C=C2)CC=3C1=CC(O)=CC3O)C4=CC(O)=C(O)C(O)=C4 '
               'NAME: Epigallocatechin 3-O-vanillate REASON: CORRECT Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'O[C@@H]1[C@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@@H]1[C@@H](O)[C@@H](O2)OC(=O)c1cc(O)c(O)c(O)c1 '
               'NAME: corilagin REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O1[C@]2([C@@](OC(=O)C=3C=4C=5C(C(OC4C(O)=C(O)C3)=O)=CC(O)=C(OC=6C(C(OC2)=O)=CC(O)=C(O)C6O)C5O)([C@@H](OC(=O)C7=CC(O)=C(O)C(O)=C7)[C@@H](OC(=O)C8=CC(O)=C(O)C(O)=C8)[C@@H]1OC(=O)C9=CC(O)=C(O)C(O)=C9)[H])[H] '
               'NAME: Bicornin REASON: CORRECT Molecule matches tannin profile '
               'with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O1[C@]2([C@@](OC(=O)C=3C(C=4C(C(OC2)=O)=CC(O)=C(O)C4O)=C(O)C(O)=C(O)C3)(C=5OC(=O)C=6C(C=7C(C(OC5C1=O)=O)=CC(O)=C(O)C7O)=C(O)C(O)=C(O)C6)[H])[H] '
               'NAME: Emblicanin B REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O[C@H]1[C@H]2[C@H](CC(O)=O)C(=O)O[C@@H]3[C@@H](COC(=O)c4cc(O)c(O)c(O)c4)O[C@@H](OC(=O)c4cc(O)c(O)c(O)c4)[C@H](OC(=O)c4cc(O)c(O)c(OC1=O)c24)[C@H]3OC(=O)c1cc(O)c(O)c(O)c1 '
               'NAME: Chebulinic acid REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](OC(=O)c4cc(O)c(O)c(O)c4)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: procyanidin B3 3-O-gallate REASON: CORRECT Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'Oc1cc(cc(Oc2c(O)c(O)c(O)cc2C(=O)O[C@H]2O[C@@H]3COC(=O)c4cc(O)c(O)c(O)c4-c4c(O)c(O)c(O)cc4C(=O)O[C@H]3[C@@H]3OC(=O)c4cc(O)c(O)c(O)c4-c4c(O)c(O)c(O)cc4C(=O)O[C@@H]23)c1O)C(=O)O[C@H]1O[C@@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@H]2[C@@H]2OC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@@H]12 '
               'NAME: agrimoniin REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'COc1cc2[C@H](OC(=O)c3ccccc3)[C@@](C)(O)[C@@H](C)Cc3cc4OCOc4c(OC)c3-c2c(OC)c1OC '
               'NAME: Schisantherin A REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H]2OC(=O)c3c1c(O)c(O)c(O)c3-c1c(O)c(O)c(O)cc1C(=O)O[C@H]2[C@@H]1OC(=O)c2cc(O)c(O)c(O)c2-c2c(O)c(O)c(O)cc2C(=O)OC[C@H]1OC(=O)c1cc(O)c(O)c(O)c1 '
               'NAME: Casuarinin REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=CC(=C1)CCCCC)C2=C(O[C@@H]3O[C@H](C(=O)O)[C@H](O)[C@@H]([C@H]3O)O)C=C(O)C=C2CCCCC '
               'NAME: Ascotricin B REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'C1[C@@H]([C@H](OC2=CC(=CC(=C21)OC(=O)C3=CC=CC=C3)OC(=O)C4=CC=CC=C4)C5=CC(=C(C=C5)OC(=O)C6=CC=CC=C6)OC(=O)C7=CC=CC=C7)OC(=O)C8=CC=CC=C8 '
               'NAME: benzoic acid '
               '[(2R,3S)-5,7-dibenzoyloxy-2-(3,4-dibenzoyloxyphenyl)-3,4-dihydro-2H-1-benzopyran-3-yl] '
               'ester REASON: CORRECT Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'S(OC1C2OC(=O)C=3C(C=4C(C(OCC1OC(O)C2O)=O)=CC(O)=C(O)C4O)=C(O)C(O)=C(O)C3)(O)(=O)=O '
               'NAME: '
               '{6,7,8,11,12,13,21,22-octahydroxy-3,16-dioxo-2,17,20-trioxatetracyclo[17.3.1.0?,?.0??,??]tricosa-4,6,8,10,12,14-hexaen-23-yl}oxidanesulfonic '
               'acid REASON: CORRECT Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=C(C(=O)O)C(=C1)C)C2=C(O)C=C(O[C@@H]3O[C@H](C(=O)O)[C@@H](O)[C@@H]([C@H]3O)O)C=C2CCCCCCCCCCCCCCC '
               'NAME: CRM646-A REASON: CORRECT Molecule matches tannin profile '
               'with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'Oc1cc(cc(O)c1O)C(=O)O[C@@H]1O[C@@H]2COC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)cc3C(=O)O[C@H]2[C@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1 '
               'NAME: eugeniin REASON: CORRECT Molecule matches tannin profile '
               'with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=CC(=C1)CCCCC)C2=C(O[C@@H]3O[C@H](C(=O)O)[C@H](O)[C@@H]([C@H]3O)O)C=C(O)C=C2CCC '
               'NAME: Ascotricin A REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O1C(OC2=C(CC=C)C=CC(OC(=O)C3=CC=CC=C3)=C2OC)C(O)C(O)C(O)C1C(O)=O '
               'NAME: '
               '6-[3-(benzoyloxy)-2-methoxy-6-(prop-2-en-1-yl)phenoxy]-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: CORRECT Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=C(C(=O)N[C@H](C(=O)O)C)C(=C1)C)C2=C(O[C@H]3O[C@H]([C@H](O)[C@H]([C@@H]3O)O)CO)C=C(OC(=O)C4=C(O)C=C(OC)C=C4C)C=C2C '
               'NAME: Amidepsine G REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'S(OC=1C(O)=CC(C(OC2C(OC(=O)C3=CC(O)=C(O)C(O)=C3)C(O)C(OC2CO)OC(=O)C4=CC(O)=C(O)C(O)=C4)=O)=CC1O)(O)(=O)=O '
               'NAME: '
               '[2,6-dihydroxy-4-({[5-hydroxy-2-(hydroxymethyl)-4,6-bis(3,4,5-trihydroxybenzoyloxy)oxan-3-yl]oxy}carbonyl)phenyl]oxidanesulfonic '
               'acid REASON: CORRECT Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC1=CC(O)=C(C(=O)N[C@H](C(=O)O)C)C(=C1)C)C2=C(O[C@@H]3O[C@@H]([C@@H](O)[C@@H]([C@H]3OC(=O)C)O)CO)C=C(OC(=O)C4=C(O)C=C(OC)C=C4C)C=C2C '
               'NAME: Amidepsine I REASON: CORRECT Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@H]3[C@H](OC(=O)c4cc(O)c(O)c(O)c4)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: procyanidin B4 3-O-gallate REASON: CORRECT Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               'False positives: SMILES: '
               'OC1=CC=C2C3(OC(C4=C3C=C(C=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O.OC1=CC=C2C3(OC(C4=C3C=CC(=C4)C(O)=O)=O)C=5C(OC2=C1)=CC(=CC5)O '
               'NAME: 5(6)-carboxyfluorescein REASON: WRONGLY CLASSIFIED '
               'Molecule matches tannin profile with multiple phenolic-like '
               'structures and appropriate linkages\n'
               ' * SMILES: '
               'C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N '
               'NAME: Cycloguanil pamoate REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'S(C1=C(C=2C(=O)C3=C(C=C(C)C=C3O)OC2C=C1O)C(=O)OC)C4=C(C=5C(=O)C6=C(C=C(C)C=C6O)OC5C=C4O)C(=O)OC '
               'NAME: Castochrin REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CC1=CNC2=C1C=CC=C2)N[C@H](C(=O)N[C@H](C(=O)O)C(C)C)CC3=CNC4=C3C=CC=C4 '
               'NAME: Trp-Trp-Val REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)CC2=C(C=CC(=C2)NC(=O)CC3=CC=CC=C3)O[C@H]1CN(C)CC4=CC=C(C=C4)C(=O)O)[C@H](C)CO '
               'NAME: '
               '4-[[[(2R,3S)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-9-[(1-oxo-2-phenylethyl)amino]-2,3,4,7-tetrahydro-1,5-benzoxazonin-2-yl]methyl-methylamino]methyl]benzoic '
               'acid REASON: WRONGLY CLASSIFIED Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'ClC1=C(OC)C=CC(=C1)C[C@H]2NC(=O)C=CC[C@H](OC(=O)[C@@H](OC([C@@H](CNC2=O)C)=O)CC(C)C)C[C@H]3O[C@@H]3C4=CC=CC=C4 '
               'NAME: Cryptophycin-40 REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'ClC1=C(O)C2=C([C@H]3C(=O)C=4C=5[C@@]2([C@H]3OC(=O)C)COC(C5C(O)=CC4C)=O)C=6C1=CC(O)=CC6C '
               'NAME: Macrosporusone C REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'CCOC(=O)C1=C(SC=C1C2=CC=C(C=C2)C)NC(=O)C3=CC4=CC=CC=C4OC3=O '
               'NAME: '
               '4-(4-methylphenyl)-2-[[oxo-(2-oxo-1-benzopyran-3-yl)methyl]amino]-3-thiophenecarboxylic '
               'acid ethyl ester REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: COC1=CC=C(C=C1)CC(=O)NC(CC(=O)O)C2=CC=C(C=C2)OC '
               'NAME: '
               '3-(4-methoxyphenyl)-3-[[2-(4-methoxyphenyl)-1-oxoethyl]amino]propanoic '
               'acid REASON: WRONGLY CLASSIFIED Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC2=CC=CC=C2)CC3=CC=CC=C3 '
               'NAME: Phe-Phe-Pro REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O[C@@H]([C@H](NC(=O)[C@@H](N)CC1=CC=CC=C1)C(=O)N[C@@H](CC=2C=3C(NC2)=CC=CC3)C(O)=O)C '
               'NAME: Phe-Thr-Trp REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'COC(=O)C1=CC=C(C=C1)COC2=CC=CC=C2C=CC(=O)C3=CC=C(C=C3)OCC(=O)O '
               'NAME: '
               '2-[4-[3-[2-[(4-methoxycarbonylphenyl)methoxy]phenyl]-1-oxoprop-2-enyl]phenoxy]acetic '
               'acid REASON: WRONGLY CLASSIFIED Molecule matches tannin '
               'profile with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: Cc1cnc2c(C(O)=O)c(Cl)ccc2c1 NAME: quinmerac REASON: '
               'WRONGLY CLASSIFIED Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: CC(=CC(=O)C)NC(CC1=CNC2=CC=CC=C21)C(=O)O NAME: '
               '3-(1H-indol-3-yl)-2-(4-oxopent-2-en-2-ylamino)propanoic acid '
               'REASON: WRONGLY CLASSIFIED Molecule matches tannin profile '
               'with multiple phenolic-like structures and appropriate '
               'linkages\n'
               ' * SMILES: '
               'O=C1O[C@@H](CC=2C1=C(O)C3=C(O)C(C4=C(O)C=5C(=O)C=C(C)C(C5C=C4OC)=O)=C(OC)C=C3C2)C '
               'NAME: Xanthoradone B REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: CCOC(=O)C1CCN(CC1)CC2=CC(=C(C=C2)OC)OCC3=CC=CC=C3 '
               'NAME: '
               '1-[(4-methoxy-3-phenylmethoxyphenyl)methyl]-4-piperidinecarboxylic '
               'acid ethyl ester REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'CN(C)S(=O)(=O)C1=CC=C(C=C1)C(=O)OCC(=O)C2=CC3=C(C=C2)OCC(=O)N3 '
               'NAME: 4-(dimethylsulfamoyl)benzoic acid '
               '[2-oxo-2-(3-oxo-4H-1,4-benzoxazin-6-yl)ethyl] ester REASON: '
               'WRONGLY CLASSIFIED Molecule matches tannin profile with '
               'multiple phenolic-like structures and appropriate linkages\n'
               ' * SMILES: '
               'O=C(OC1=C(C(O)=C(C(=O)O)C(=C1)C)C)C2=C(O)C(=C(OC)C=C2C)C=O '
               'NAME: Baeomycesic acid REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'C12=C(C=CC=C1)NC=3[C@]([C@@]4(/C(/CN(CC4)CCC32)=C\\C)[H])(C(=O)OC)COC(C)=O '
               'NAME: O-acetyl-15alpha-stemmadenine REASON: WRONGLY CLASSIFIED '
               'Molecule matches tannin profile with multiple phenolic-like '
               'structures and appropriate linkages\n'
               ' * SMILES: '
               'O=C(N[C@@H](CC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCC(O)=O)CC=1C=2C(NC1)=CC=CC2 '
               'NAME: Glu-Trp-Asp REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=[N+]([O-])C1=C(C(OC2=C(C(O)=CC(=C2)C)C(=O)OC[C@@H](C(=O)OC)C)=C(OC)C=C1O)C(=O)OC '
               'NAME: Pseudogymnoascin A REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'C(=O)(C=1N(C)C(=CC1)CC(NCC(OC2=C(OC)C=CC=C2)=O)=O)C3=CC=C(C=C3)C '
               'NAME: amtolmetin guacil REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               ' * SMILES: '
               'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)CCC(=O)N '
               'NAME: Trp-Gln-Gln REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC3=CC=CC=C3)[C@H](CC)C '
               'NAME: Phe-Ile-Trp REASON: WRONGLY CLASSIFIED Molecule matches '
               'tannin profile with multiple phenolic-like structures and '
               'appropriate linkages\n'
               ' * SMILES: '
               'C1=2C=3C=4[C@](N(C[C@@H](C4)C(N[C@@H](CC)CO)=O)C)(CC1=CN(C2C=CC3)C)[H].C(=C\\C(O)=O)\\C(O)=O '
               'NAME: Methysergide maleate REASON: WRONGLY CLASSIFIED Molecule '
               'matches tannin profile with multiple phenolic-like structures '
               'and appropriate linkages\n'
               'False negatives: SMILES: '
               'O=C1OC(=CC=2C1=C(O)C(OC)=C(O)C2C3=C(O)C(OC)=C(O)C4=C3C=C(C)OC4=O)C '
               'NAME: Bireticulol REASON: MISSED Lacking typical tannin '
               'linkages\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1cc(O)c(O)c(O)c1 '
               'NAME: Gallocatechin-(4alpha->8)-epigallocatechin REASON: '
               'MISSED Lacking typical tannin linkages\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1cc(O)c(O)c(O)c1 '
               'NAME: (+)-catechin-(4alpha->8)-(-)-epigallocatechin REASON: '
               'MISSED Lacking typical tannin linkages\n'
               ' * SMILES: COc1cc2CC(C)C(C)Cc3cc(OC)c(OC)c(OC)c3-c2c(OC)c1OC '
               'NAME: Deoxyschizandrin REASON: MISSED Lacking typical tannin '
               'linkages\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1OC(=O)C=2C=C(OC)C(OC)=C(O)C2)CO '
               'NAME: [3,4,5-Trihydroxy-6-(hydroxymethyl)oxan-2-yl] '
               '3-hydroxy-4,5-dimethoxybenzoate REASON: MISSED Found only 1 '
               'phenolic-like structures, need at least 2\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)ccc34)c3cc(O)c(O)c(O)c3)c(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)ccc34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: '
               'Robinetinidol-(4alpha->8)-catechin-(6->4alpha)-robinetinidol '
               'REASON: MISSED Lacking typical tannin linkages\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1cc(O)c(O)c(O)c1 '
               'NAME: (-)-epicatechin-(4beta->8)-(+)-gallocatechin REASON: '
               'MISSED Lacking typical tannin linkages\n'
               ' * SMILES: O1C2=C3C=4C(=CC(OC)=C(O)C4OC(=O)C3=CC(O)=C2OC)C1=O '
               'NAME: 2,8-Di-O-methylellagic acid REASON: MISSED Lacking '
               'typical tannin linkages\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4c3c(O)cc3O[C@@]5(Oc6cc(O)cc(O)c6[C@@H]([C@H]5O)c43)c3ccc(O)c(O)c3)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: cinnamtannin B-1 REASON: MISSED Lacking typical tannin '
               'linkages\n'
               ' * SMILES: '
               'COc1cc2C[C@@H](C)[C@@H](C)[C@@H](O)c3cc4OCOc4c(OC)c3-c2c(OC)c1OC '
               'NAME: Isogomisin O REASON: MISSED Lacking typical tannin '
               'linkages\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Found only 1 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Found only 1 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Found only 1 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Found only 1 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Found only 0 aromatic rings, '
                                               'need at least 2'}],
    'sample_false_negatives': [   {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(=O)C=2C=C(OC)C(OC)=C(O)C2)CO',
                                      'name': '[3,4,5-Trihydroxy-6-(hydroxymethyl)oxan-2-yl] '
                                              '3-hydroxy-4,5-dimethoxybenzoate',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'COC1=CC(=CC=C1O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C(O)=O)O)O)C)C(O)=O',
                                      'name': 'Vanillic acid-4-O-glucuronide',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(OC(=O)C=2C=C(OC3OC(C(O)C(O)C3O)C(O)=O)C(O)=C(O)C2)C(O)C1O)CO',
                                      'name': '6-[2,3-dihydroxy-5-({[2,3,5-trihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy}carbonyl)phenoxy]-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC=2C(O)=C(C=C(C2)C(OC)=O)CC=C(C)C)CO',
                                      'name': 'Methyl '
                                              '3,4-dihydroxy-5-prenylbenzoate '
                                              '3-glucoside',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(=O)C2=CC(O)=C(O)C(O)=C2)CO',
                                      'name': '[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl] '
                                              '3,4,5-trihydroxybenzoate',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(=O)C2=CC(OC)=C(OC)C(OC)=C2)CO',
                                      'name': '1-O-(3,4,5-Trimethoxybenzoyl)-b-D-glucopyranoside',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC2=C(OC)C=C(C=C2)C(O)=O)CO',
                                      'name': 'Vanillic acid glucoside',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC=2C(O)=C(C(O)=CC2O)C(O)=O)C(O)=O',
                                      'name': '6-(3-carboxy-2,4,6-trihydroxyphenoxy)-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O[C@@H]1C[C@](O)(C[C@@H](OC(=O)\\C=C\\c2ccc(O)c(O)c2)[C@@H]1O)C(O)=O',
                                      'name': 'chlorogenic acid',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'},
                                  {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(=O)C2=CC(OC)=C(O)C=C2)C(O)=O',
                                      'name': '3,4,5-trihydroxy-6-(4-hydroxy-3-methoxybenzoyloxy)oxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Found only 1 aromatic rings, '
                                                'need at least 2'}],
    'attempt': 3,
    'reasoning': 'To improve the classification of tannin molecules using '
                 'their SMILES strings, I will focus on refining the detection '
                 'of phenolic-like structures and the complex linkage patterns '
                 "typical of tannins. I'll incorporate more detailed "
                 'substructure searches, including specific aromatic rings and '
                 'linkage patterns such as catechin units and ellagic '
                 'derivatives. \n'
                 '\n'
                 'Here is a step-by-step plan before implementing the code:\n'
                 '\n'
                 '1. **Phenolic and Aromatic Detection**: Ensure detection of '
                 'multiple aromatic or phenolic rings. I need to refine the '
                 'SMARTS pattern to detect potential tannin-related phenolic '
                 'structures more accurately.\n'
                 '\n'
                 '2. **Linkage Patterns**: Enhanced recognition of relevant '
                 'ester and ether bonds that are typical in tannins. Consider '
                 'more specific linker patterns that reflect complex tannins, '
                 'including bonds between different polyphenolic units.\n'
                 '\n'
                 '3. **Intra- and Intermolecular Links**: Consider special '
                 'linkage combinations between phenolic units that facilitate '
                 'branching, common in hydrolyzable and condensed tannins.\n'
                 '\n'
                 '4. **Implementation in RDKit**: Use RDKit for updated SMARTS '
                 'patterns and confirm presence of specific linkages and ring '
                 'systems more accurately.\n'
                 '\n'
                 '5. **Validation with Logical Flow**: Ensure all logic steps '
                 'aggregate to connect recognized structures into verified '
                 'tannin classification.\n'
                 '\n'
                 "Let's proceed with the new implementation by making some "
                 'strategic adjustments to the substructure search to better '
                 'capture tannin characteristics:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 144,
    'num_false_positives': 100,
    'num_true_negatives': 1076,
    'num_false_negatives': 25,
    'num_negatives': None,
    'precision': 0.5901639344262295,
    'recall': 0.8520710059171598,
    'f1': 0.6973365617433414,
    'accuracy': 0.9070631970260223,
    'negative_predictive_value': 0.9772933696639419}