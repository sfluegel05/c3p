"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by either:
    1. A benzene ring with oxygen-containing functional groups and long hydrocarbon chains
    2. A fatty acid/ethanolamide structure with long hydrocarbon chains

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - cannabinoids typically have higher molecular weights
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a cannabinoid"

    # Define patterns for cannabinoid features
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ethanolamide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    cannabinoid_pattern = Chem.MolFromSmarts("c1cc(O)c2[C@@H]3C=C(CC[C@H]3C(C)(C)Oc2c1)")

    # Check for classic cannabinoid structure (benzene + oxygen + long chain)
    has_benzene = mol.HasSubstructMatch(benzene_pattern)
    has_oxygen = len(mol.GetSubstructMatches(oxygen_pattern)) >= 1
    has_long_chain = mol.HasSubstructMatch(long_chain_pattern)

    # Check for ethanolamide/fatty acid structure
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_fatty_acid = mol.HasSubstructMatch(fatty_acid_pattern)

    # Check for specific cannabinoid pattern
    has_cannabinoid_pattern = mol.HasSubstructMatch(cannabinoid_pattern)

    # Classification logic
    if has_cannabinoid_pattern:
        return True, "Contains specific cannabinoid structural motif"
    elif (has_benzene and has_oxygen and has_long_chain):
        return True, "Contains a benzene ring with oxygen-containing functional groups and long hydrocarbon chains"
    elif (has_ethanolamide or has_fatty_acid) and has_long_chain:
        return True, "Contains ethanolamide/fatty acid structure with long hydrocarbon chain"
    elif has_benzene and has_oxygen:
        return False, "Contains benzene and oxygen but lacks characteristic long chain"
    else:
        return False, "Does not match cannabinoid structural patterns"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:67194',
                          'name': 'cannabinoid',
                          'definition': 'A diverse group of pharmacologically '
                                        'active secondary metabolite '
                                        'characteristic to Cannabis plant as '
                                        'well as produced naturally in the '
                                        'body by humans and animals. '
                                        'Cannabinoids contain oxygen as a part '
                                        'of the heterocyclic ring or in the '
                                        'form of various functional groups. '
                                        'They are subdivided on the basis of '
                                        'their origin.',
                          'parents': ['CHEBI:36963'],
                          'xrefs': [   'PMID:22234284',
                                       'PMID:22530636',
                                       'PMID:22555283',
                                       'Wikipedia:Cannabinoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.3067484662576688 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCCCN1C=C(C(=O)CC2=CC=C(OC)C=C2)C2=CC=CC=C12 NAME: JWH-201 '
               'REASON: CORRECT Contains a benzene ring with oxygen-containing '
               'functional groups and long hydrocarbon chains\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C '
               'NAME: Roselipin 3D REASON: CORRECT Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: O(C(CO)CO)C(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC '
               'NAME: MG(0:0/20:3(5Z,8Z,11Z)/0:0) REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'CCCCCc1cc(O)c2[C@@H]3C=C(CC[C@H]3C(C)(C)Oc2c1)C(O)=O NAME: '
               '11-nor-9-carboxy-Delta(9)-tetrahydrocannabinol REASON: CORRECT '
               'Contains a benzene ring with oxygen-containing functional '
               'groups and long hydrocarbon chains\n'
               ' * SMILES: '
               '[H][C@@]12CCC(C)=C[C@@]1([H])c1c(O)c(C(O)=O)c(CCCCC)cc1OC2(C)C '
               'NAME: Delta(9)-tetrahydrocannabinolic acid REASON: CORRECT '
               'Contains a benzene ring with oxygen-containing functional '
               'groups and long hydrocarbon chains\n'
               ' * SMILES: '
               '[H][C@]1(CCC(C)=C[C@H]1C1=C(O)C=C(CCCCC)C=C1O)C(C)=C NAME: '
               'cannabidiol REASON: CORRECT Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'CCCCCCC(C)(C)C=1C=CC([C@@]2(C[C@](CC[C@]2(CCCO)[H])([H])O)[H])=C(C1)O '
               'NAME: CP-55940 REASON: CORRECT Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O[C@@H](CCCCC)/C=C/C=C\\C/C=C\\C/C=C\\CCCC(OC(CO)CO)=O NAME: '
               '15-HETE-G REASON: CORRECT Contains ethanolamide/fatty acid '
               'structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@H](O)[C@H](O)CO)[C@H](O)CO)/C(=C/[C@@H]([C@H](O)/C(=C/[C@@H]([C@H](O)/C(=C/[C@@H]([C@H](O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)/C(=C/[C@H](C[C@H](CC)C)C)/C)C)/C)C)/C)C)/C '
               'NAME: TMC-171A REASON: CORRECT Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               '[H]C(CO)(CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC NAME: '
               '2-arachidonoylglycerol REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'O(C(C(O)([2H])[2H])(C(O)([2H])[2H])[2H])C(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC '
               'NAME: 2-arachidonoyl glycerol-d5 REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: CCCCCC=CCC=CCC=CCC=CCCCC(=O)OC(CO)CO NAME: '
               'eicosa-5,8,11,14-tetraenoic acid 1,3-dihydroxypropan-2-yl '
               'ester REASON: CORRECT Contains ethanolamide/fatty acid '
               'structure with long hydrocarbon chain\n'
               ' * SMILES: CCCCCc1cc(O)c(C\\C=C(\\C)CCC=C(C)C)c(O)c1C(O)=O '
               'NAME: cannabinerolic acid REASON: CORRECT Contains a benzene '
               'ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C(CO)CO NAME: '
               'MG(0:0/20:2(11Z,14Z)/0:0) REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\CCCCCCCC)C(CO)CO NAME: '
               'MG(0:0/20:1(11Z)/0:0) REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               '[H][C@@]12CCC(C)=C[C@@]1([H])c1c(O)cc(CCCCC)cc1OC2(C)C NAME: '
               'Delta(9)-tetrahydrocannabinol REASON: CORRECT Contains a '
               'benzene ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'O=C(O[C@@H]([C@@H](O)[C@H](O)CO)CO)C(=CC(C(O)C(=CC(C(O)C(=CC(C(O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@@H]1O)O)CO)C(=CC(CC(CC)C)C)C)C)C)C)C)C)C '
               'NAME: Roselipin 3E REASON: CORRECT Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: CCCCCc1cc(O)c(C\\C=C(/C)CCC=C(C)C)c(O)c1C(O)=O '
               'NAME: cannabigerolic acid REASON: CORRECT Contains a benzene '
               'ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: CCCCCc1cc(O)c2[C@@H]3C=C(CO)CC[C@H]3C(C)(C)Oc2c1 '
               'NAME: 11-hydroxy-Delta(9)-tetrahydrocannabinol REASON: CORRECT '
               'Contains a benzene ring with oxygen-containing functional '
               'groups and long hydrocarbon chains\n'
               ' * SMILES: O(C(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:3(11Z,14Z,17Z)/0:0) REASON: CORRECT Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               '[H][C@]1(C=C(C)CC[C@H]1C(C)=C)c1c(O)cc(CCCCC)c(C(O)=O)c1O '
               'NAME: cannabidiolic acid REASON: CORRECT Contains a benzene '
               'ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C(CO)CO '
               'NAME: MG(0:0/20:4(8Z,11Z,14Z,17Z)/0:0) REASON: CORRECT '
               'Contains ethanolamide/fatty acid structure with long '
               'hydrocarbon chain\n'
               ' * SMILES: '
               '[C@@H]1(C=C([C@@H]2C([C@H]1C2)(C)C)CO)C=3C(=CC(=CC3OC)C(CCCCCC)(C)C)OC '
               'NAME: HU-308 REASON: CORRECT Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O(C(CO)CO)C(=O)CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC '
               'NAME: 12-HETE-G REASON: CORRECT Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: C1(=C(C=C(CCCCC)C=C1O)O)C/C=C(/CCC=C(C)C)\\C NAME: '
               'cannabigerol REASON: CORRECT Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               'False positives: SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(O)(O)=O '
               'NAME: PA(18:3(6Z,9Z,12Z)/20:1(11Z)) REASON: WRONGLY CLASSIFIED '
               'Contains ethanolamide/fatty acid structure with long '
               'hydrocarbon chain\n'
               ' * SMILES: '
               'C1CC(C1)CNC[C@H]2[C@@H]([C@@H](N2C(=O)CC3=CN=CC=C3)CO)C4=CC=CC=C4 '
               'NAME: '
               '1-[(2R,3S,4R)-2-[(cyclobutylmethylamino)methyl]-4-(hydroxymethyl)-3-phenyl-1-azetidinyl]-2-(3-pyridinyl)ethanone '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O=C(O[C@H](CC(=O)N)CCCCCCCCCCC(C)C)[C@@H](NC(=O)C[C@@H](O)CCCCCCCCCCCC(C)C)CCCC(=O)O '
               'NAME: WA8242A1 REASON: WRONGLY CLASSIFIED Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'COC1=CC=CC(=C1)C(=O)N2C[C@@H]3[C@@H]([C@H](N3CC4=CC=CC=C4Cl)CO)C5=CC=CC=C52 '
               'NAME: '
               '[(1S,2aS,8bS)-2-[(2-chlorophenyl)methyl]-1-(hydroxymethyl)-1,2a,3,8b-tetrahydroazeto[2,3-c]quinolin-4-yl]-(3-methoxyphenyl)methanone '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O=C1C2=C(C(=O)C(=C1C(C)C)C)C3[C@@H](C4=CC=CC=C4)C56C2(C(=O)C(=C(O)C5=O)C(C)C)C6[C@H]3C7=CC=CC=C7 '
               'NAME: Carbocyclinone-534 REASON: WRONGLY CLASSIFIED Contains a '
               'benzene ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CC(C)C)C(O)=O)CO '
               'NAME: Cys-Ser-Leu REASON: WRONGLY CLASSIFIED Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               '[H][C@]12CC[C@]([H])([C@H]([C@H](C1)OC(=O)c1ccccc1)C(=O)OC)[NH+]2C '
               'NAME: cocaine(1+) REASON: WRONGLY CLASSIFIED Contains a '
               'benzene ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'C=12[C@]3(C)[C@@]([H])(CCC1C=C(C(=C2)O)C(C)C)C(C[C@@H](C3)O)(C)C '
               'NAME: salviol REASON: WRONGLY CLASSIFIED Contains a benzene '
               'ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'C[C@@H]1CCCCO[C@H]([C@H](CN(C(=O)C2=C(O1)C=CC(=C2)NS(=O)(=O)C)[C@H](C)CO)C)CN(C)C '
               'NAME: '
               'N-[(3R,9R,10S)-9-[(dimethylamino)methyl]-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-13-oxo-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-16-yl]methanesulfonamide '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O=C(N[C@@H](CCC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)CC=1C=2C(NC1)=CC=CC2 '
               'NAME: Asp-Trp-Glu REASON: WRONGLY CLASSIFIED Contains a '
               'benzene ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'P(O[C@H]1[C@H](O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)(OCC(OC(=O)C)COC(=O)C)(O)=O '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-2-{[(1R,2R,3R,4R,5R,6S)-2-({[2,3-bis(acetyloxy)propoxy](hydroxy)phosphoryl}oxy)-3,4,5,6-tetrahydroxycyclohexyl]oxy}-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]ethanimidic '
               'acid REASON: WRONGLY CLASSIFIED Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'P(OCC(OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCCCC)(OCCNC)(O)=O '
               'NAME: '
               '{2-[(7Z,10Z,13Z,16Z)-docosa-7,10,13,16-tetraenoyloxy]-3-(tetracosanoyloxy)propoxy}[2-(methylamino)ethoxy]phosphinic '
               'acid REASON: WRONGLY CLASSIFIED Contains ethanolamide/fatty '
               'acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'C1C[C@H]([C@H](O[C@@H]1CCNC(=O)C2=CC3=C(C=C2)OCO3)CO)NS(=O)(=O)C4=CC=C(C=C4)Cl '
               'NAME: '
               'N-[2-[(2S,5R,6S)-5-[(4-chlorophenyl)sulfonylamino]-6-(hydroxymethyl)oxan-2-yl]ethyl]-1,3-benzodioxole-5-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O=C1OC([C@@H]2CCC=3[C@]4([C@]([C@@H]([C@@H]([C@H]5OC(=O)C(C)=CC5)C)C[C@H]4OC(=O)C)(C)CC(C3C=C2C=C1)=O)C)(C)C '
               'NAME: 11-oxo-colobetaolactone E REASON: WRONGLY CLASSIFIED '
               'Contains ethanolamide/fatty acid structure with long '
               'hydrocarbon chain\n'
               ' * SMILES: O=C1O[C@@H](CCCCCCCCCCCCC)C(=C1C)C(=O)O NAME: '
               'Lichesterinic acid REASON: WRONGLY CLASSIFIED Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'COC(=O)C[C@@H]1C[C@H]2[C@@H]([C@@H](O1)CO)OC3=C2C=C(C=C3)NS(=O)(=O)C '
               'NAME: methyl '
               '2-[(1S,3S,4aR,9aS)-1-(hydroxymethyl)-6-(methanesulfonamido)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b][1]benzofuran-3-yl]acetate '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C=2C1=C(O)C3=C(C(=O)C4(O)C5C3OC4(C(O[C@@H]6O[C@@H]([C@@H](O[C@H]7O[C@@H]([C@@](O)(C(=O)C)[C@H](C7)O)C)CC6)C)C(C5OC8O[C@H]([C@H](O)[C@@H](C8)O)C)O)CCC)C2O)C=9C(OC)=C(OC)C%10=CC=%11CC%12C(O[C@@H]%13O[C@H]([C@H](O)[C@@H](C%13)O)C)C(O)C(C(C%12(C(C%11C(=C%10C9O)O)=O)O)(O)CCC)O[C@@H]%14O[C@@H]([C@@H](O[C@H]%15O[C@@H]([C@@](O)(C(=O)C)[C@H](C%15)O)C)CC%14)C '
               'NAME: Hibarimicin G REASON: WRONGLY CLASSIFIED Contains a '
               'benzene ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=CC=CC(=C2OCC[C@@H]3CC[C@@H]([C@@H](O3)CN(C1=O)C)OC)C#N)C '
               'NAME: LSM-9211 REASON: WRONGLY CLASSIFIED Contains a benzene '
               'ring with oxygen-containing functional groups and long '
               'hydrocarbon chains\n'
               ' * SMILES: '
               'O1C([C@@H](O)C(O)C(O)[C@@H]1OC2=C(OC=3C(C2=O)=C(O)C=C(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)C3)C5=CC(OC)=C(O)C=C5)CO[C@@H]6OC([C@H](O)[C@H](O)C6O)C '
               'NAME: isorhamnetin 3-rutinoside-7-glucoside REASON: WRONGLY '
               'CLASSIFIED Contains a benzene ring with oxygen-containing '
               'functional groups and long hydrocarbon chains\n'
               ' * SMILES: O(CCCCCCCCCCC)C(=O)C(N(C)C(OC)=O)C NAME: Alanine, '
               'N-methyl-N-methoxycarbonyl-, undecyl ester REASON: WRONGLY '
               'CLASSIFIED Contains ethanolamide/fatty acid structure with '
               'long hydrocarbon chain\n'
               ' * SMILES: '
               'O=C1NC=2C(=O)C(=C(OC)C(C2)=O)[C@H](O)[C@H](C[C@H](OC)[C@H](O)[C@H](C=C([C@@H]([C@H](C=CC=C1C)OC)OC(=O)N)C)C)C '
               'NAME: 15-Hydroxygeldanamycin REASON: WRONGLY CLASSIFIED '
               'Contains ethanolamide/fatty acid structure with long '
               'hydrocarbon chain\n'
               ' * SMILES: '
               'CN(C)C(=O)C[C@@H]1C[C@@H]2[C@H]([C@H](O1)CO)OC3=C2C=C(C=C3)NS(=O)(=O)C4=CC=CC(=C4)OC '
               'NAME: '
               '2-[(1R,3S,4aS,9aR)-1-(hydroxymethyl)-6-[(3-methoxyphenyl)sulfonylamino]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N,N-dimethylacetamide '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O '
               'NAME: '
               '(4aS,6S,7S,7aR)-6-[6-amino-8-[(4-chlorophenyl)thio]-9-purinyl]-2-hydroxy-2-oxo-4a,6,7,7a-tetrahydro-4H-furo[3,2-d][1,3,2]dioxaphosphorin-7-ol '
               'REASON: WRONGLY CLASSIFIED Contains a benzene ring with '
               'oxygen-containing functional groups and long hydrocarbon '
               'chains\n'
               ' * SMILES: '
               'ClC(C(O)CCCCCCCC)C/C=C/C=C\\[C@H](O)CC(=O)NC=1C(=O)O[C@](C1)(CC)C '
               'NAME: Korormicin H REASON: WRONGLY CLASSIFIED Contains '
               'ethanolamide/fatty acid structure with long hydrocarbon chain\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OC(CCCCCCCCCCC)CC(=O)N[C@H]1[C@@H](OP(O)(O)=O)O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H](OC(=O)CC(CCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)[C@H]2NC(=O)CC(CCCCCCCCCCC)OC(=O)CCCCCCCCCCC)[C@@H](O)[C@@H]1OC(=O)CC(O)CCCCCCCCCCC '
               'NAME: diphospho heptaacyl lipid A REASON: WRONGLY CLASSIFIED '
               'Contains ethanolamide/fatty acid structure with long '
               'hydrocarbon chain\n'
               'False negatives: SMILES: CCCCCCCCCCCCCCCC(=O)NCCO NAME: '
               'palmitoyl ethanolamide REASON: MISSED Does not match '
               'cannabinoid structural patterns\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCCl NAME: '
               "arachidonyl-2'-chloroethylamide REASON: MISSED Does not match "
               'cannabinoid structural patterns\n'
               ' * SMILES: C(CCC)C/C=C\\C/C=C\\CC1C(C/C=C\\CCCC(NCCO)=O)O1 '
               'NAME: N-[(5Z,11Z,14Z)-8,9-epoxyicosatrienoyl]ethanolamine '
               'REASON: MISSED Does not match cannabinoid structural patterns\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)NCCO NAME: '
               'anandamide REASON: MISSED Does not match cannabinoid '
               'structural patterns\n'
               ' * SMILES: N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF '
               'NAME: 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole REASON: '
               'MISSED Does not match cannabinoid structural patterns\n'
               ' * SMILES: C(CCCCCCCCCCCCC/C=C\\CCCCCCCC)(NCCO)=O NAME: '
               'N-(15Z-tetracosenoyl)ethanolamine REASON: MISSED Does not '
               'match cannabinoid structural patterns\n'
               ' * SMILES: '
               'Cc1c(C(=O)c2cccc3ccccc23)c2cccc3OC[C@@H](CN4CCOCC4)n1c23 NAME: '
               'WIN 55212-2 REASON: MISSED Contains benzene and oxygen but '
               'lacks characteristic long chain\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCOC(CO)CO NAME: '
               '2-arachidonyl glyceryl ether REASON: MISSED Does not match '
               'cannabinoid structural patterns\n'
               ' * SMILES: C(CCCO)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(NCCO)=O '
               'NAME: N-(20-hydroxyarachidonoyl)ethanolamine REASON: MISSED '
               'Does not match cannabinoid structural patterns\n'
               ' * SMILES: CCCCCCCC\\C=C/CCCCCCCC(=O)NCCO NAME: oleoyl '
               'ethanolamide REASON: MISSED Does not match cannabinoid '
               'structural patterns\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cl.C=1C=C(OC)C(OC)=C(C1CN2CCNCC2)OC.Cl',
                                     'name': 'Trimetazidine hydrochloride',
                                     'reason': 'Contains benzene and oxygen '
                                               'but lacks characteristic long '
                                               'chain'},
                                 {   'smiles': '[H]S(=O)(=O)C[C@H](N)C(O)=O',
                                     'name': 'L-cysteine-S-dioxide',
                                     'reason': 'Molecular weight too low for a '
                                               'cannabinoid'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C=C(CO[C@@H]1CN(C)CC3=CC=CC=C3OC)N=N2)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'Contains benzene and oxygen '
                                               'but lacks characteristic long '
                                               'chain'},
                                 {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'Does not match cannabinoid '
                                               'structural patterns'},
                                 {   'smiles': 'O=C1[C@H]([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@H](O)CC(=C)C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)C',
                                     'name': 'Gilvsin C',
                                     'reason': 'Does not match cannabinoid '
                                               'structural patterns'},
                                 {   'smiles': 'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C(=O)C',
                                     'name': '1-[(4S,7S,8S)-5-acetyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-phenylurea',
                                     'reason': 'Contains benzene and oxygen '
                                               'but lacks characteristic long '
                                               'chain'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3=NC=CN=C3',
                                     'name': 'N-[(4S,7R,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[oxo(2-pyrazinyl)methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Contains benzene and oxygen '
                                               'but lacks characteristic long '
                                               'chain'},
                                 {   'smiles': 'COc1ccc(cc1)\\C=C(NC=O)\\C(NC=O)=C\\c1ccc(OS(O)(=O)=O)cc1',
                                     'name': 'Fumiformamide',
                                     'reason': 'Contains benzene and oxygen '
                                               'but lacks characteristic long '
                                               'chain'},
                                 {   'smiles': 'CCCCC(CC)COCCCN',
                                     'name': '3-(2-ethylhexoxy)propan-1-amine',
                                     'reason': 'Molecular weight too low for a '
                                               'cannabinoid'},
                                 {   'smiles': 'N(=CC1=CC=CC=C1)CCC2=CC=CC=C2',
                                     'name': 'N-benzylidene-2-phenylethanamine',
                                     'reason': 'Does not match cannabinoid '
                                               'structural patterns'}],
    'sample_false_negatives': [   {   'smiles': 'N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF',
                                      'name': '1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole',
                                      'reason': 'Does not match cannabinoid '
                                                'structural patterns'},
                                  {   'smiles': 'Cc1c(C(=O)c2cccc3ccccc23)c2cccc3OC[C@@H](CN4CCOCC4)n1c23',
                                      'name': 'WIN 55212-2',
                                      'reason': 'Contains benzene and oxygen '
                                                'but lacks characteristic long '
                                                'chain'},
                                  {   'smiles': 'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCOC(CO)CO',
                                      'name': '2-arachidonyl glyceryl ether',
                                      'reason': 'Does not match cannabinoid '
                                                'structural patterns'}],
    'attempt': 2,
    'reasoning': 'The previous code had several issues that led to both false '
                 'positives and false negatives:\n'
                 '\n'
                 '1. **False Positives**: The code incorrectly classified '
                 'molecules that had benzene rings with oxygen-containing '
                 'functional groups and long hydrocarbon chains, even if they '
                 "weren't cannabinoids. This is because the criteria were too "
                 "broad and didn't account for the specific structural motifs "
                 'found in cannabinoids.\n'
                 '\n'
                 '2. **False Negatives**: The code missed several true '
                 'cannabinoids, particularly those with ethanolamide/fatty '
                 'acid structures, because the patterns used to identify these '
                 'structures were not specific enough. For example, the '
                 'ethanolamide pattern `[NX3][CX3](=[OX1])CCO` was too '
                 'restrictive and missed many valid cannabinoids.\n'
                 '\n'
                 '3. **Molecular Weight Filter**: The molecular weight filter '
                 'was too restrictive, potentially excluding smaller '
                 'cannabinoids.\n'
                 '\n'
                 '### Improvements:\n'
                 '1. **Refine Patterns**: Use more specific patterns to '
                 'identify cannabinoid structures. For example, the '
                 'ethanolamide pattern should be more inclusive to capture a '
                 'wider range of cannabinoids.\n'
                 '\n'
                 '2. **Combine Patterns**: Use combined patterns to ensure '
                 'that the molecule has both the necessary functional groups '
                 'and the correct structural arrangement.\n'
                 '\n'
                 '3. **Adjust Molecular Weight Filter**: Lower the molecular '
                 'weight threshold to include smaller cannabinoids.\n'
                 '\n'
                 '4. **Add Additional Checks**: Include checks for specific '
                 'functional groups and structural motifs that are '
                 'characteristic of cannabinoids.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 35,
    'num_false_positives': 100,
    'num_true_negatives': 132,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.25925925925925924,
    'recall': 0.9210526315789473,
    'f1': 0.40462427745664736,
    'accuracy': 0.6185185185185185,
    'negative_predictive_value': 0.9777777777777777}