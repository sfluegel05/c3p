"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Isoflavonoids
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Expand isoflavonoid definition: 1-benzopyran core with diverse aryl group
    isoflavonoid_patterns = [
        Chem.MolFromSmarts("Oc1ccc2ccccc2c1-c3ccccc3"),         # Typical isoflavonoid
        Chem.MolFromSmarts("O=c1cc2ccccc2oc1-c3ccccc3"),        # Variant with additional ketone
        Chem.MolFromSmarts("Oc1ccc2ccccc2c1-c3cc(O)ccc3"),      # With hydroxyl groups
        Chem.MolFromSmarts("Oc1ccc2cc(O)cccc2c1-c3ccc(OC)cc3")  # With methoxy groups
    ]

    for pattern in isoflavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches isoflavonoid substructure with possible functional variations"
    
    # Check possible higher complexity with glucose or other sugar attachments
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)O)O)CO")
    if mol.HasSubstructMatch(sugar_pattern):
        return True, "Contains potential isoflavonoid structure with sugar attachments"

    # Check for common substituents in isoflavonoids
    common_substituents = [
        Chem.MolFromSmarts("CO"),  # Methoxy group
        Chem.MolFromSmarts("C=C"), # Carbon double bonds
        Chem.MolFromSmarts("CC(=O)O") # Acetate or ester linkages
    ]

    for substituent in common_substituents:
        if mol.HasSubstructMatch(substituent):
            return True, "Contains functional group common in isoflavonoids"

    return False, "Does not contain the necessary isoflavonoid structure or functional groups"

# Example validation
examples = [
    "COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O",  # 2,3-dihydrobiochanin A
    "CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O"  # neobavaisoflavone
]

for e in examples:
    result, reason = is_isoflavonoid(e)
    print(f"SMILES: {e}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50753',
                          'name': 'isoflavonoid',
                          'definition': 'Any 1-benzopyran with an aryl '
                                        'substituent at position 3. The term '
                                        'was originally restricted to natural '
                                        'products, but is now also used to '
                                        'describe semi-synthetic and fully '
                                        'synthetic compounds.',
                          'parents': [   'CHEBI:26004',
                                         'CHEBI:38443',
                                         'CHEBI:72544'],
                          'xrefs': [   'MetaCyc:Isoflavonoids',
                                       'Wikipedia:Isoflavonoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.8172484599589322 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: COc1ccc(cc1)C1COc2cc(O)cc(O)c2C1=O '
               'NAME: 2,3-dihydrobiochanin A REASON: CORRECT Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O(C1C(O)C(O)C(OC1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(OC)C=C4)CO)C5OCC(O)(C5O)CO '
               'NAME: Formononetin 7-O-apiosyl-(1->2)-glucoside REASON: '
               'CORRECT Contains isoflavonoid structure with methoxy group\n'
               ' * SMILES: [H][C@]1(COc2cc(C)ccc2[C@H]1O)c1ccccc1OC NAME: '
               'conferol B REASON: CORRECT Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=C(OC)C=C3C(OC=C(C3=O)C4=CC=C(OC)C=C4)=C2)CO '
               'NAME: Wistin REASON: CORRECT Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'Oc1ccc-2c(c1)[C@@H]1Oc3ccccc3[C@H]3CC(=O)c4c(O)cc(O)c-2c4[C@@H]13 '
               'NAME: Ohioensin-A REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'O1C2=C(C(=O)C(C3=CC=C(OC)C=C3)=C1)C=CC(OCC=C(C)C)=C2 NAME: '
               'Maximaisoflavone J REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: O=C1[C@H](COc2ccccc12)c1ccccc1 NAME: '
               '(3S)-isoflavanone REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'O1C(C=CC2=C1C=CC(C=3C(=O)C=4C(O)=C(CCC(O)(C)C)C(O)=CC4OC3)=C2O)(C)C '
               'NAME: Kanzonol T REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: COc1cc(O)c2c(c1)occ(-c1ccc(O)cc1)c2=O NAME: '
               'prunetin REASON: CORRECT Contains isoflavonoid structure with '
               'methoxy group\n'
               ' * SMILES: COc1ccc(cc1)-c1coc2cc(OC)c(OC)c(OC)c2c1=O NAME: '
               "5,6,7,4'-tetramethoxyisoflavone REASON: CORRECT Contains "
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: Oc1ccc(cc1)-c1coc2cc3OCOc3c(O)c2c1=O NAME: irilone '
               'REASON: CORRECT Contains isoflavonoid structure with methoxy '
               'group\n'
               ' * SMILES: '
               'O1C(C=CC2=C1C=C(O)C3=C2OC=C(C4=CC5=C(OC(C=C5)(C)C)C=C4)C3=O)(C)C '
               'NAME: Ulexone B REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: O=C1C2=C(OC[C@@]1(O)C3=CC=C(O)C=C3)C=C(O)C=C2O '
               'NAME: Dothideoisoflavanone REASON: CORRECT Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=CC=3OCC(CC3C=C2)C4=C(O)C(OC)=C(OC)C=C4)CO '
               'NAME: Isomucronulator 7-O-glucoside REASON: CORRECT Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O=C1C(C2=CC=C(O[C@@H]3O[C@H]([C@H](O)[C@H]([C@H]3O)OC)C)C=C2)=COC4=C1C=CC(=C4)O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)OC)C '
               'NAME: Termisoflavone C REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'O1C(C(O)CC=2C1=C(C=CC2O)C=3C(=O)C4=C(OC3)C=C(O)C(=C4O)CC=C(C)C)(C)C '
               'NAME: Lupinisolone C REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=3OCC(CC3C=C2)C4=C(O)C=C(OC)C=C4)CO '
               'NAME: Vestitone 7-glucoside REASON: CORRECT Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O1[C@@H](OC2=CC=3OC[C@@H](CC3C=C2)C4=CC=C(O)C=C4)[C@H](O)[C@@H](O)[C@H](O)[C@H]1C(O)=O '
               'NAME: Equol 7-O-glucuronide REASON: CORRECT Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1C(O)=O)C(C(OCC2OC(OC=3C=C4OC=C(C(=O)C4=C(O)C3)C5=CC=C(OC)C=C5)C(O)C(O)C2O)=O)C(OC)=O '
               'NAME: '
               '3,4,5-trihydroxy-6-{1-methoxy-1,3-dioxo-3-[(3,4,5-trihydroxy-6-{[5-hydroxy-3-(4-methoxyphenyl)-4-oxo-4H-chromen-7-yl]oxy}oxan-2-yl)methoxy]propan-2-yl}oxane-2-carboxylic '
               'acid REASON: CORRECT Contains isoflavonoid structure with '
               'methoxy group\n'
               ' * SMILES: O1C(C(C2=C1C=C3OCC(CC3=C2)C4=C(O)C=C(O)C=C4)(C)C)C '
               'NAME: Cyclomillinol REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: O=C1C2=C(OC=C1C3=CC=C(O)C=C3)C(OC)=C(O)C=C2O NAME: '
               'Psi-tectorigenin REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: C1Oc2ccccc2C[C@H]1c1ccccc1 NAME: (S)-isoflavan '
               'REASON: CORRECT Contains isoflavonoid structure with methoxy '
               'group\n'
               ' * SMILES: '
               'O1C(C=CC=2C1=C(C=3OC=C(C(=O)C3C2O)C4=CC(O)=C(O)C=C4)CC=C(C)C)(C)C '
               'NAME: Auriculasin REASON: CORRECT Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'CC(C)=CCc1c(O)c2c(occ(-c3ccc(O)cc3)c2=O)c2C=CC(C)(C)Oc12 NAME: '
               'Osajin REASON: CORRECT Contains isoflavonoid structure with '
               'methoxy group\n'
               ' * SMILES: '
               'O1CC(CC2=C1C3=C(OC(C=C3)(C)C)C=C2)C4=C(OC)C=C(O)C=C4 NAME: '
               "2'-O-Methylglabridin REASON: CORRECT Contains isoflavonoid "
               'structure with methoxy group\n'
               'False positives: SMILES: '
               'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2 '
               'NAME: Besigomsin REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1 '
               'NAME: (+)-tephrorin B REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O '
               'NAME: ancitabine hydrochloride REASON: WRONGLY CLASSIFIED '
               'Contains isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=CC(=C3)C#N)O[C@H]1CN(C)C(=O)C4=CC=NC=C4)[C@H](C)CO '
               'NAME: '
               'N-[[(4R,5R)-8-(3-cyanophenyl)-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-N-methyl-4-pyridinecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)NC(=O)NC(C)C)[C@H](C)CO)C)CN(C)C '
               'NAME: '
               '1-[(3R,9S,10R)-9-[(dimethylamino)methyl]-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-16-yl]-3-propan-2-ylurea '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc3O[C@@]4(C[C@@H](c5c(O)cc(O)cc5O4)c3c2O[C@@H]1c1ccc(O)cc1)c1ccc(O)cc1 '
               'NAME: mahuannin D REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'CN1[C@@H]([C@H]2CCN([C@H]2C3=C1C=CC(=C3)C#CC4=CC=C(C=C4)OC)C(=O)NC5=CC=C(C=C5)C(F)(F)F)CO '
               'NAME: '
               '(3aS,4S,9bR)-4-(hydroxymethyl)-8-[2-(4-methoxyphenyl)ethynyl]-5-methyl-N-[4-(trifluoromethyl)phenyl]-3,3a,4,9b-tetrahydro-2H-pyrrolo[3,2-c]quinoline-1-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C3=CC=C(C=C3)OC)C(=O)N(C[C@@H]1OC)C)C)C '
               'NAME: '
               '4-methoxy-N-[(4S,7S,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]benzamide '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'C[C@@H]1CCCCO[C@H]([C@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@H](C)CO)C)CN(C)CC3=CC=C(C=C3)C(=O)O '
               'NAME: '
               '4-[[[(3R,9R,10S)-16-(dimethylamino)-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-9-yl]methyl-methylamino]methyl]benzoic '
               'acid REASON: WRONGLY CLASSIFIED Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'C1=C(C=C2C(=C1OC)O[C@@H]([C@H]2C(=O)[O-])C3=CC=C(C(=C3)OC)O)/C=C/C([O-])=O '
               'NAME: (+)-DCA-CC(2-) REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@@H]3[C@@H](CC[C@@H](O3)CCN(C1=O)C)OC)C '
               'NAME: LSM-38104 REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'CN1C[C@@H]2[C@@H](CC[C@@H](O2)CCOC3=C(C=C(C=C3)C#N)C(=O)N(C[C@H](CC1=O)C4=CC=CC=C4)C)OC '
               'NAME: LSM-38314 REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: O=C1[C@]2(C(C3C(CC2)C=4C(CC3)=CC(OC)=C(O)C4)CC1)C '
               'NAME: 2-Hydroxy-3-methoxyestrone REASON: WRONGLY CLASSIFIED '
               'Contains isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC(C)C)C(=O)N(C[C@@H]1OC)C)C)CC3=CC=C(C=C3)C4=CC=CC=N4 '
               'NAME: '
               '1-[(4S,7S,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[[4-(2-pyridinyl)phenyl]methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-propan-2-ylurea '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@@H]1CN(C)S(=O)(=O)C4=CC=CC=C4OC)[C@H](C)CO '
               'NAME: LSM-8955 REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O=C1C(=C(NC(=C1)OC)C/C=C(/C/C=C/C(=C/C(C(O)/C(=C/C)/C)C)/C)\\C)C '
               'NAME: Mer-A2026B REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'O(C1CC(O)(CC2=C1C(O)=C3C(=C2O)C([O-])C4=C(C3[O-])C(OC)=CC=C4)C(=O)CO)C5OC(C(O)C(N)C5)C '
               'NAME: Doxorubicin-semiquinone REASON: WRONGLY CLASSIFIED '
               'Contains isoflavonoid structure with methoxy group\n'
               ' * SMILES: O1C=2C(=CC(=C(OC)C2)C(O)=O)C=CC1=O NAME: Buntansin '
               'A REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: '
               'COC1=CC=CC=C1C(=O)OCC#CCSC2=NN=C(O2)C3=CC=CC4=CC=CC=C43 NAME: '
               '2-methoxybenzoic acid '
               '4-[[5-(1-naphthalenyl)-1,3,4-oxadiazol-2-yl]thio]but-2-ynyl '
               'ester REASON: WRONGLY CLASSIFIED Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'COC1=CC=CC=C1S(=O)(=O)N2C[C@H](COC[C@H]3[C@H]2CC[C@H](O3)CC(=O)OC)O '
               'NAME: '
               '2-[(3R,6aR,8S,10aR)-3-hydroxy-1-(2-methoxyphenyl)sulfonyl-3,4,6,6a,8,9,10,10a-octahydro-2H-pyrano[2,3-c][1,5]oxazocin-8-yl]acetic '
               'acid methyl ester REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'CC=CC1=CN=C2C(=C1)C(=O)N(C[C@@H]([C@@H](O2)CN(C)C(=O)NC3=CC=C(C=C3)F)C)[C@@H](C)CO '
               'NAME: LSM-7398 REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: CN1CCc2cc3OCOc3c3c2[C@@H]1CC31C=CC(=O)C=C1 NAME: '
               'Mecambrine REASON: WRONGLY CLASSIFIED Contains isoflavonoid '
               'structure with methoxy group\n'
               ' * SMILES: '
               'O1C2=C(O)C(O)=CC(=C2C[C@H]3[C@]1(CC=CC(C[C@@H]4[C@](CC3)(OC5=C(C(=CC(=C5)O)C)C4)C)(C)C)C)C '
               'NAME: 22-hydroxyramiferin REASON: WRONGLY CLASSIFIED Contains '
               'isoflavonoid structure with methoxy group\n'
               ' * SMILES: '
               'CCCC#CC1=CC2=C(N=C1)O[C@@H]([C@@H](CN(C2=O)[C@H](C)CO)C)CN(C)C(=O)NC3=CC(=CC=C3)F '
               'NAME: '
               '3-(3-fluorophenyl)-1-[[(2S,3R)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-8-pent-1-ynyl-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-2-yl]methyl]-1-methylurea '
               'REASON: WRONGLY CLASSIFIED Contains isoflavonoid structure '
               'with methoxy group\n'
               ' * SMILES: CC(=O)OC(C)(C)C1CC2=C(O1)C=C3C(=C2)C=CC(=O)O3 NAME: '
               'acetic acid '
               '2-(7-oxo-2,3-dihydrofuro[3,2-g][1]benzopyran-2-yl)propan-2-yl '
               'ester REASON: WRONGLY CLASSIFIED Contains isoflavonoid '
               'structure with methoxy group\n'
               'False negatives: SMILES: '
               'ClC1=C2OC=C(C3=CC=C(O)C=C3)C(C2=C(Cl)C(=C1O)O)=O NAME: '
               'Maduraktermol G REASON: MISSED Does not contain the necessary '
               'isoflavonoid structure or functional groups\n'
               ' * SMILES: CC(C)=CCc1cc(ccc1O)-c1coc2cc(O)ccc2c1=O NAME: '
               'neobavaisoflavone REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'O=C1C=2C=3C=4C(OC=5C=CC=C(C5C4[C@@H](O)[C@@H](C3C(=C1)O)C)O)=CC2O '
               'NAME: (-)-dalmanol A REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C(C(O)=C3C2=C4C(=C5C(=O)C=CC=C5OC4=C1)[C@H]([C@@H](CC(=O)C6=C(O)C=CC=C6O)C)[C@H]3C)=C(O)C '
               'NAME: (+)-acetodalmanol B REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'O1C2=C(C/C=C(/CO)\\C)C(O)=CC(O)=C2C(=O)C(C3=C(O)C=C(O)C=C3)=C1 '
               'NAME: 2,3-Dehydrokievitol REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'C=1(C2C(C=CC1O)C(C(=CO2)C=3C=CC(=CC3)O)=O)[C@@H]4O[C@@H]([C@H](C([C@H]4O)=O)O)CO '
               "NAME: 3''-oxopuearin REASON: MISSED Does not contain the "
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'S(OC1=CC=C(C=2C(=O)C3=C(OC2)C=C(O)C=C3)C=C1)(O)(=O)=O NAME: '
               "Daidzein 4'-sulfate REASON: MISSED Does not contain the "
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'O1C2=C(C3=C1OC4=C(C3=O)C(O)=C(C(O)=C4)CC=C(C)C)C=CC(O)=C2CC=C(C)C '
               'NAME: Lupinalbin F REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               ' * SMILES: '
               'O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=C(O)C=C(O)C=C3)=C1 '
               "NAME: 2',4',5,7-Tetrahydroxy-8-prenylisoflavone REASON: MISSED "
               'Does not contain the necessary isoflavonoid structure or '
               'functional groups\n'
               ' * SMILES: '
               'O1C2=C(C/C=C(/CO)\\C)C(O)=CC(O)=C2C(=O)C(=C1)C3=CC=C(O)C=C3 '
               'NAME: Gancaonin C REASON: MISSED Does not contain the '
               'necessary isoflavonoid structure or functional groups\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'OC1=C(NC(=O)NC2=CC=CC=C2Br)C=CC(=C1)[N+]([O-])=O',
                                     'name': 'SB225002',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'C1CNC(=O)C2=C1C3=C(C=CC(=C3N2)Cl)Cl',
                                     'name': '5,8-dichloro-2,3,4,9-tetrahydropyrido[3,4-b]indol-1-one',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'Cl[Si](Cl)(Cl)CCl',
                                     'name': 'Trichloro(chloromethyl)silane',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'CCCCN1CCCC[C@H]1C(=O)Nc1c(C)cccc1C',
                                     'name': 'levobupivacaine',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'CC(CCCCCCCC=O)C',
                                     'name': '9-methyldecanal',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'O=C1OC(=CC(=C1C)O)CCCC(C)C',
                                     'name': 'Violapyrone A',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)NC2CC2',
                                     'name': 'N-cyclopropyl-4-methylbenzenesulfonamide',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'},
                                 {   'smiles': 'O=C(C=1NC(=O)NC1CC)C2=CC=NC=C2',
                                     'name': 'piroximone',
                                     'reason': 'Does not contain the necessary '
                                               'isoflavonoid structure or '
                                               'functional groups'}],
    'sample_false_negatives': [   {   'smiles': 'ClC1=C2OC=C(C3=CC=C(O)C=C3)C(C2=C(Cl)C(=C1O)O)=O',
                                      'name': 'Maduraktermol G',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'S(OC1=CC=C(C=2C(=O)C3=C(OC2)C=C(O)C=C3)C=C1)(O)(=O)=O',
                                      'name': "Daidzein 4'-sulfate",
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'O1C2=C(C(=O)C(C3=CC=C(O)C=C3)=C1)C=CC(O)=C2O',
                                      'name': '8-Hydroxydaidzein',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'S(OC1=CC=C(C=2C(=O)C=3C(OC2)=CC(O)=CC3O)C=C1)(O)(=O)=O',
                                      'name': '[4-(5,7-dihydroxy-4-oxo-4H-chromen-3-yl)phenyl]oxidanesulfonic '
                                              'acid',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'O1C2=C(C(=O)C(C3=CC(O)=C(O)C=C3)=C1)C(O)=CC(O)=C2O',
                                      'name': "5,7,8,3',4'-Pentahydroxyisoflavone",
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'Oc1cc(O)c2c(c1)occ(-c1ccc(O)c(O)c1)c2=O',
                                      'name': 'orobol',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'BrC1=C2OC=C(C3=CC(Br)=C(O)C=C3)C(C2=C(O)C(=C1O)Br)=O',
                                      'name': 'Maduraktermol N',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'O1C=2C(C(=O)C(C3=CC=CC=C3)=C1)=C(O)C=C(O)C2',
                                      'name': '5,7-Dihydroxyisoflavone',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'O=C1C=2C=3C=4C(OC=5C=CC=C(C5C4[C@H]([C@@H](CC(=O)C6=C(O)C=CC=C6O)C)[C@H](C3C(=C1)O)C)O)=CC2O',
                                      'name': '(+)-dalmanol B',
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'},
                                  {   'smiles': 'O=C1C(C2=CC(O)=CC(=C2)O)=COC3=C1C=CC(=C3)O',
                                      'name': "7-hydroxy-3',5'-dihydroxyisoflavone",
                                      'reason': 'Does not contain the '
                                                'necessary isoflavonoid '
                                                'structure or functional '
                                                'groups'}],
    'attempt': 3,
    'reasoning': 'To improve the classification of isoflavonoids using the '
                 "SMILES string, let's first analyze the shortcomings of the "
                 'previous implementation. The `is_isoflavonoid` function '
                 'primarily checked for a specific pattern of the isoflavonoid '
                 'class - the 1-benzopyran core with an aryl group at position '
                 '3. However, this approach failed to consider the broader '
                 'diversity of isoflavonoid structures, ending up with a mix '
                 'of false positives and false negatives.\n'
                 '\n'
                 'One key issue is that isoflavonoids come in various forms '
                 'with additional substituents not limited to methoxy or '
                 'hydroxyl groups, leading to the misclassification of many '
                 'true positives and false positives. Therefore, we should:\n'
                 '\n'
                 '1. Adjust the substructure search to better capture the '
                 'benzopyran core and its connectivity with the 3-position '
                 'aryl group.\n'
                 '2. Handle potential isoflavonoid variations, such as '
                 'different functional groups or synthetic modifications.\n'
                 '\n'
                 "Here's how we can update the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 443,
    'num_false_positives': 100,
    'num_true_negatives': 18,
    'num_false_negatives': 33,
    'num_negatives': None,
    'precision': 0.8158379373848987,
    'recall': 0.930672268907563,
    'f1': 0.8694798822374878,
    'accuracy': 0.7760942760942761,
    'negative_predictive_value': 0.35294117647058826}