"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound composed of two or more pyrrole units.
    Pyrrole units are five-membered rings containing exactly one nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()

    # Initialize count of pyrrole units
    pyrrole_units = 0

    # Iterate over rings
    for ring_atoms in ring_info.AtomRings():
        # Check if ring is size 5
        if len(ring_atoms) != 5:
            continue

        # Count nitrogen atoms in the ring
        num_nitrogen = 0
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                num_nitrogen += 1

        # If ring has exactly one nitrogen, count it as a pyrrole unit
        if num_nitrogen == 1:
            pyrrole_units += 1

    # Check if molecule is a polypyrrole
    if pyrrole_units >= 2:
        return True, f"Contains {pyrrole_units} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_units} pyrrole unit(s). A polypyrrole requires two or more pyrrole units"


__metadata__ = {
    'chemical_class': {
        'name': 'polypyrrole',
        'definition': 'A compound composed of two or more pyrrole units.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38077',
                          'name': 'polypyrrole',
                          'definition': 'A compound composed of two or more '
                                        'pyrrole units.',
                          'parents': ['CHEBI:38101'],
                          'xrefs': ['Beilstein:8538310'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7807017543859649 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'Cc1c(CCC(O)=O)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(C)c5C=C)c(C)c4C=C)c(C)c3CCC(O)=O '
               'NAME: protoporphyrin REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'CC1=C(CCC(O)=O)C2=[N+]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N+]5[Mg--]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2 '
               'NAME: magnesium protoporphyrin REASON: CORRECT Contains 2 '
               'pyrrole units\n'
               ' * SMILES: '
               'CC1=C2CC3=C(C(=C(N3)C=C4C5=CC=C(C([C@]5(C(=CC6=NC(=CC(=C1CCC(=O)O)N2)C(=C6C)CCC(=O)O)N4)C)C(=O)OC)C(=O)OC)C)C=C '
               'NAME: LSM-4252 REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               '[H]C(C)=C1C(C)C2=CC3=C(C(C)=O)C(C)=C4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(N34)N3C(=CC1=[N+]25)C(C)=C1C(=O)[C@H](C(=O)OC)C6=C31 '
               'NAME: bacteriochlorophyll b REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'C1=2N3C(C=C4N5=C(C=C6N7C(=CC8=N(C(=C1)C(=C8CCC(O)=O)C([H])=O)[Fe]735)C(=C6C)CCC(O)=O)C(=C4C)C=C)=C(C2C)[C@H](CC/C=C(\\C)/CC/C=C(/CCC=C(C)C)\\C)O '
               'NAME: ferroheme a REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'CC1(C)c2ccc([nH]2)C(C)(C)c2ccc([nH]2)C(C)(C)c2ccc([nH]2)C(C)(C)c2ccc1[nH]2 '
               'NAME: meso-octamethylporphyrinogen REASON: CORRECT Contains 4 '
               'pyrrole units\n'
               ' * SMILES: '
               'O=C1C=2N=C(C1(C3OC(C(O)CC3O)C)C)C=C4NC(=CC5=NC(=CC=6NC(C2)=C(C6)C)C(C5=O)(C7OC(C(O)CC7OC(=O)C)C)C)C=C4C '
               'NAME: Tolyporphin C REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'c1ccc2c3nc(nc4[nH]c(nc5[nH]c(n3)c3ccccc53)c3ccccc43)c2c1 NAME: '
               'subphthalocyanine REASON: CORRECT Contains 3 pyrrole units\n'
               ' * SMILES: '
               'Bc1c2ccc(n2)c(B)c2ccc([nH]2)c(B)c2ccc(n2)c(B)c2ccc1[nH]2 NAME: '
               '5,10,15,20-tetraborylporphyrin REASON: CORRECT Contains 2 '
               'pyrrole units\n'
               ' * SMILES: '
               'O(C=1C(N=C(C1)C=2NC=CC2)=CC=3NC=4CCCCCCCC[C@@H](C3C4)CC(C)C)C '
               'NAME: Prodigiosin R1 REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               '[H]C(=C([H])C1=[N+]2C(C=C1)=Cc1ccc(-c3ccc[nH]3)n1[B-]2(F)F)c1ccc(OCC(=O)NCCCCCC(=O)ON2C(=O)CCC2=O)cc1 '
               'NAME: BODIPY 650/665-X REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               '[H]C(=O)c1c(CC)c2cc3[nH]c4c([C@@H](C(=O)OC)C(=O)c4c3C)c3nc(cc4[nH]c(cc1n2)c(C=C)c4C)[C@@H](C)[C@@H]3CCC(O)=O '
               'NAME: pheophorbide b REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'OC(=O)CCC=1C=2NC(C1C)=CC3=NC(C(=C3CCC(O)=O)C)=CC=4NC(=C(C4C=C)C)C=C5N=C(C(=C5C)CCC(O)=O)C2 '
               'NAME: Harderoporphyrin REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'COC(=O)CCc1c(C)c2=CC3=[N]4C(=Cc5c(C=C)c(C)c6C=C7C(C)=C(CCC(O)=O)C8=[N]7[Mg]4(n56)n2c1=C8)C(C)=C3C=C '
               'NAME: magnesium protoporphyrin 13-monomethyl ester REASON: '
               'CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'CCc1c(C)c2[nH]c1cc1nc(cc3[nH]c(cc4nc(cc5[nH]c2c(C)c5CC)c(C)c4CCCOC(=O)N(CCOCCOCCOC)CCOCCOCCOC)c(CC)c3CC)c(CCCOC(=O)N(CCOCCOCCOC)CCOCCOCCOC)c1C '
               'NAME: sapphyrin PCI-2050 REASON: CORRECT Contains 3 pyrrole '
               'units\n'
               ' * SMILES: Bc1c2ccc(cc3ccc(cc4ccc(cc5ccc1[nH]5)n4)[nH]3)n2 '
               'NAME: 5-borylporphyrin REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'CCC1=C(CC)C2=[N]3C1=Cc1c(CC)c(CC)c4C=C5[N]6=C(C=c7c(CC)c(CC)c(=C2)n7[Pt]36n14)C(=O)C5(CC)CC '
               'NAME: platinum(II) octaethylporphyrin ketone REASON: CORRECT '
               'Contains 2 pyrrole units\n'
               ' * SMILES: '
               '[H]C(=O)c1c(C)c2\\C=C3/N=C([C@@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)[C@@H]3C)C3=c4c(C(=O)[C@@H]3C(=O)OC)c(C)c3=CC5=N\\C(=C/c1n2[Mg]n43)\\C(C)=C5CC '
               'NAME: chlorophyll d REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: Brc1c2ccc(cc3ccc(cc4ccc(cc5ccc1[nH]5)n4)[nH]3)n2 '
               'NAME: 10-bromoporphyrin REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'C=1(\\C=C/2\\NC(=O)C(=C2C=C)C)NC(=C(C1C)CCC(=O)O)CC=3NC(\\C=C/4\\NC(C(=C4C)C=C)=O)=C(C3CCC(=O)O)C '
               'NAME: (E,E)-bilirubin REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'CCC1=C(CO)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31 '
               'NAME: 7(1)-hydroxychlorophyll a REASON: CORRECT Contains 2 '
               'pyrrole units\n'
               ' * SMILES: '
               'CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(O)=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31 '
               'NAME: chlorophyllide a REASON: CORRECT Contains 2 pyrrole '
               'units\n'
               ' * SMILES: '
               'CCC1=C(C)C(=O)NC1Cc1[nH]c(Cc2[nH]c(CC3NC(=O)C(C=C)=C3C)c(C)c2CCC(O)=O)c(CCC(O)=O)c1C '
               'NAME: D-Urobilinogen REASON: CORRECT Contains 2 pyrrole units\n'
               ' * SMILES: '
               'CC[C@@H]1[C@@H](C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(O)=O)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31 '
               'NAME: 3-vinylbacteriochlorophyllide a REASON: CORRECT Contains '
               '2 pyrrole units\n'
               ' * SMILES: '
               'CCC1=C(C=O)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31 '
               'NAME: chlorophyll b REASON: CORRECT Contains 2 pyrrole units\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'C1=CC(=CC2=NC(C=C2)=CC2=NC(C=C2)=Cc2ccc[nH]2)N=C1 NAME: '
               '21H-bilin REASON: MISSED Contains only 1 pyrrole unit(s). A '
               'polypyrrole requires two or more pyrrole units\n'
               ' * SMILES: '
               'OC(=O)CCC1C(C=2NC1=C(C3=NC(C(=C3C(O)=O)C)=CC=4N=C(C=C5N=C(C2)C(=C5C=C)C)/C(/C4CC)=C\\O)CC(O)=O)C '
               'NAME: Rhodin G7 REASON: MISSED Contains only 0 pyrrole '
               'unit(s). A polypyrrole requires two or more pyrrole units\n'
               ' * SMILES: '
               'C[C@@]1(CCC(O)=O)[C@@H](CC(O)=O)C2=C3N4C(=CC5=[N+]6C(C[C@]7(C)N8C(=CC1=[N+]2[Co--]468)C(CCC(O)=O)=C7CC(O)=O)=C(CCC(O)=O)[C@]5(C)CC(O)=O)[C@@H](CCC(O)=O)[C@]3(C)CC(O)=O '
               'NAME: cobalt-precorrin-5B REASON: MISSED Contains only 0 '
               'pyrrole unit(s). A polypyrrole requires two or more pyrrole '
               'units\n'
               ' * SMILES: '
               'C[C@]1(CC(O)=O)[C@H](CCC(O)=O)C2=CC3=[N+]4C(=Cc5c(CCC(O)=O)c(CC(O)=O)c6C=C7[N+]8=C(C=C1N2[Co--]48n56)[C@@H](CCC(O)=O)[C@]7(C)CC(O)=O)C(CCC(O)=O)=C3CC(O)=O '
               'NAME: cobalt-sirohydrochlorin REASON: MISSED Contains only 1 '
               'pyrrole unit(s). A polypyrrole requires two or more pyrrole '
               'units\n'
               ' * SMILES: '
               'Cc1ccc2C=C3C=CC(CCCCCCCCCCCC(O)=O)=[N+]3[B-](F)(F)n12 NAME: '
               'BODIPY-3823 REASON: MISSED Contains only 1 pyrrole unit(s). A '
               'polypyrrole requires two or more pyrrole units\n'
               ' * SMILES: '
               '[H][C@]12N=C(C(C)C3=N[C@@](C)(CC4=C(CCC(O)=O)[C@](C)(CC(O)=O)C(/C=C5\\N[C@]1(C)[C@@](C)(CC(O)=O)[C@@H]5CCC(O)=O)=N4)C(C)=C3CCC(O)=O)[C@](C)(CCC(O)=O)[C@H]2CC(O)=O '
               'NAME: precorrin-7 REASON: MISSED Contains only 0 pyrrole '
               'unit(s). A polypyrrole requires two or more pyrrole units\n'
               ' * SMILES: '
               'CC1OC(=O)C[C@@]2(C)[C@H](CCC(O)=O)C3=CC4=[N+]5C(=CC6=[N+]7C(=CC8=[N+]9C(=C(CC(O)=O)[C@@]8(C)CCC(O)=O)C12N3[Co--]579)C(CCC(O)=O)=C6CC(O)=O)[C@@H](CCC(O)=O)[C@]4(C)CC(O)=O '
               'NAME: cobalt(II)-factor IV REASON: MISSED Contains only 0 '
               'pyrrole unit(s). A polypyrrole requires two or more pyrrole '
               'units\n'
               ' * SMILES: '
               'C1(NC(/C=C/2\\N=C(/C=C/3\\N\\C(=C/C4=NC(C(=C4CCC(O)=O)C)=O)\\C(=C3C)CCC(O)=O)C(=C2C)C=C)=C(C1C)C=C)=O '
               'NAME: biliverdin delta REASON: MISSED Contains only 1 pyrrole '
               'unit(s). A polypyrrole requires two or more pyrrole units\n'
               ' * SMILES: '
               '[H][C@]12N=C(\\C(C)=C3/N=C(/C=C4\\N=C(\\C(C)=C5/N[C@]1(C)[C@@](C)(CC(O)=O)[C@@H]5CCC(O)=O)[C@@](C)(CC(O)=O)[C@@H]4CCC(O)=O)C(C)(C)[C@@H]3CCC(O)=O)[C@](C)(CCC(O)=O)[C@H]2CC(O)=O '
               'NAME: hydrogenobyrinic acid REASON: MISSED Contains only 0 '
               'pyrrole unit(s). A polypyrrole requires two or more pyrrole '
               'units\n'
               ' * SMILES: '
               '[H][C@]12N=C(C(C)C3=N[C@@](C)(CC4=N\\C(=C(C)/C5=N[C@]1(C)[C@@](C)(CC(O)=O)[C@@H]5CCC(O)=O)[C@@](C)(CC(O)=O)[C@@H]4CCC(O)=O)C(C)=C3CCC(O)=O)[C@](C)(CCC(O)=O)[C@H]2CC(O)=O '
               'NAME: precorrin-8X REASON: MISSED Contains only 0 pyrrole '
               'unit(s). A polypyrrole requires two or more pyrrole units\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Contains only 0 pyrrole '
                                               'unit(s). A polypyrrole '
                                               'requires two or more pyrrole '
                                               'units'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'The previous code only counted pyrrole units that are '
                 'aromatic five-membered rings with exactly one nitrogen atom '
                 'where all atoms are marked as aromatic. However, in many of '
                 'the false negatives, the pyrrole units are not aromatic '
                 "according to RDKit's perception, possibly due to differences "
                 'in tautomeric forms or resonance structures.\n'
                 '\n'
                 'To improve the classification, we can modify the code to '
                 'count all five-membered rings containing exactly one '
                 'nitrogen atom, regardless of aromaticity. This approach will '
                 'capture both aromatic and non-aromatic pyrrole units.\n'
                 '\n'
                 'Here is the revised code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 283,
    'num_false_positives': 100,
    'num_true_negatives': 3341,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.7389033942558747,
    'recall': 1.0,
    'f1': 0.8498498498498499,
    'accuracy': 0.9731471535982814,
    'negative_predictive_value': 1.0}