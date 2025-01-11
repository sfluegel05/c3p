"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signaling molecule derived from 
    the oxygenation of C20 fatty acids, excluding classic leukotrienes and prostanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure 20 carbon count (specific count for C20 fatty acid backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Does not have at least 20 carbons, which is required for C20 fatty acids"

    # Check for at least three double bonds, enforced for icosatrienoid backbones
    dbl_bond_count = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if dbl_bond_count < 3:
        return False, "Insufficient double bonds (at least three required for icosatrienes)"

    # Presence and count of essential functional groups (e.g., OH, epoxy, COOH)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4][OX2H]")))
    epoxy_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2R]1[CX3][CX3]1")))
    carboxylic_acid_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O")))
    
    if hydroxyl_count < 1 and epoxy_count < 1 and carboxylic_acid_count < 1:
        return False, "Lacking essential functional groups (need at least one OH, epoxy, or COOH group)"

    # Exclude classic leukotrienes or prostanoids
    leukotriene_exclusion = Chem.MolFromSmarts("[CH2X4]=[CHX3][CH2X4][CHX4]=[CHX3][CHX4]C") 
    prostanoid_exclusion = Chem.MolFromSmarts("C1=CCCCC1[C](=O)O")

    if mol.HasSubstructMatch(leukotriene_exclusion) or mol.HasSubstructMatch(prostanoid_exclusion):
        return False, "Matches exclusion patterns for classic icosanoids (leukotriene or prostanoid)"

    return True, "Molecule is a nonclassic icosanoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61703',
                          'name': 'nonclassic icosanoid',
                          'definition': 'Any biologically active signalling '
                                        'molecule made by oxygenation of C20 '
                                        'fatty acids other than the classic '
                                        'icosanoids (the leukotrienes and the '
                                        'prostanoids).',
                          'parents': ['CHEBI:23899'],
                          'xrefs': ['PMID:10477044'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.3178807947019867 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCCC[C@@H](O)[C@H](O)\\C=C\\C=C\\C=C/C=C/[C@@H](O)CCCC(O)=O '
               'NAME: 15-epi-lipoxin B4 REASON: CORRECT Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: O=C(CCC/C=C\\C/C=C\\CC(/C=C/C(C(CCCCC)O)O)O)O NAME: '
               '11,14,15-trihydroxy-(5Z,8Z,12E)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: C(CCCO)CC1C(C/C=C\\C/C=C\\C/C=C\\CCCC(O)=O)O1 NAME: '
               '14,15-epoxy-20-hydroxy-(5Z,8Z,11Z)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(/C=C\\C=C\\[C@H](CCCCCO)O)=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)O '
               'NAME: 20-hydroxylipoxin A4 REASON: CORRECT Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\[C@H](CCCC(O)=O)O)=C\\C=C\\C=C\\[C@@H](C\\C=C/C=C/[C@@H](CC)O)O '
               'NAME: resolvin E1 REASON: CORRECT Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: O=C(CCC/C=C\\C/C=C\\C/C=C\\C(C1C(CCCCC)O1)O)O NAME: '
               '13-hydroxy-14,15-epoxy-(5Z,8Z,11Z)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C[C@H]1O[C@@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O '
               'NAME: (8R)-hydroxy-(11R,12R)-epoxyicosa-(5Z,9E,14Z)-trienoic '
               'acid REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O '
               'NAME: (13S)-hydroxy-(14S,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic '
               'acid REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: O=C(CCC/C=C\\C/C=C\\CC(C(/C=C/C(CCCCC)O)O)O)O NAME: '
               '11,12,15-trihydroxy-(5Z,8Z,13E)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\CC(/C=C/[C@@H]1[C@H](CCCCC)O1)O)O NAME: 11 '
               'hydroxy-(14R,15S)-epoxy-(5Z,8Z,12E)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/[C@H](O)[C@H]1O[C@H]1C\\C=C/CCCC(O)=O NAME: '
               '(8S,9S)-epoxy-(10R)-hydroxyicosa-(5Z,11Z,14Z)-trienoic acid '
               'REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\C=C/C=C/C=C/[C@H](O)[C@@H](O)CCCC(O)=O '
               'NAME: 6-epi-lipoxin A4 REASON: CORRECT Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\C=C/C=C/[C@@H](C/C=C\\CC)O)=C/C=C/C(C(CCCC(=O)O)O)O NAME: '
               '(7E,9E,11Z,13E,15R,17Z)-5,6,15-trihydroxyicosapentaenoic acid '
               'REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\C=C/C=C/C=C/[C@H]([C@H](CCCCCO)O)O)=C/[C@H](CCCC(O)=O)O '
               'NAME: 20-hydroxylipoxin B4 REASON: CORRECT Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: C(CCCO)C/C=C\\C/C=C\\C/C=C\\CC1C(CCCC(O)=O)O1 NAME: '
               '5,6-epoxy-20-hydroxy-(8Z,11Z,14Z)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: C(CCCO)C/C=C\\CC1C(C/C=C\\C/C=C\\CCCC(O)=O)O1 NAME: '
               '11,12-epoxy-20-hydroxy-(5Z,8Z,14Z)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC[C@@H]1O[C@H]1[C@H](O)\\C=C/C\\C=C/C\\C=C/CCCC(O)=O NAME: '
               '(13R)-hydroxy-(14S,15S)-epoxyicosa-(5Z,8Z,11Z)-trienoic acid '
               'REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\[C@H](CCCC(O)=O)O)=C\\C=C\\C=C\\[C@@H](C/C=C\\C=C\\C(CC)=O)O '
               'NAME: 18-oxoresolvin E1 REASON: CORRECT Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\C([C@@H]1[C@H](CCCCC)O1)O)O NAME: '
               '13-hydroxy-(14R,15S)-epoxy-(5Z,8Z,11Z)-icosatrienoic acid '
               'REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: C(CCCO)C/C=C\\C/C=C\\CC1C(C/C=C\\CCCC(O)=O)O1 NAME: '
               '8,9-epoxy-20-hydroxy-(5Z,11Z,14Z)-icosatrienoic acid REASON: '
               'CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC[C@H](O)[C@H](O)\\C=C\\C=C\\C=C/C=C/[C@@H](O)CCCC(O)=O '
               'NAME: lipoxin B4 REASON: CORRECT Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: O=C(CCC/C=C\\C/C=C\\C/C=C\\C(C(C(CCCCC)O)O)O)O '
               'NAME: 13,14,15-trihydroxy-(5Z,8Z,11Z)-icosatrienoic acid '
               'REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O '
               'NAME: (8R)-hydroxy-(11S,12S)-epoxyicosa-(5Z,9E,14Z)-trienoic '
               'acid REASON: CORRECT Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'CCCCCC(=O)\\C=C\\C=C/C=C/C=C/[C@@H](O)[C@@H](O)CCCC(O)=O NAME: '
               '15-oxolipoxin A4 REASON: CORRECT Molecule is a nonclassic '
               'icosanoid\n'
               'False positives: SMILES: '
               'C(\\CCCC(O)=O)=C\\C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](CC)O NAME: '
               '18(S)-HEPE REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C=C/C(O)C\\C=C/CCCC(O)=O NAME: '
               '8-HETE REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: '
               'O(O)[C@@H](CCCC([O-])=O)/C=C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: 5(S)-HpEPE(1-) REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: C(CCC)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC([O-])=O)O '
               'NAME: 17-HETE(1-) REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'C(CCC)C/C=C\\C/C=C\\C[C@@H]1[C@H](C/C=C\\CCCC([O-])=O)O1 NAME: '
               '(8S,9R)-EET(1-) REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: C([C@@H](CC)O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(O)=O '
               'NAME: 18(R)-HETE REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: CCC(O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O NAME: '
               '18-HETE REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: '
               'C([C@H](/C=C/C=C/C=C/[C@H](CCCC(O)=O)O)O)/C=C\\CCCCCO NAME: '
               '20-hydroxy-6-trans-leukotriene B4 REASON: WRONGLY CLASSIFIED '
               'Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\C([C@H]1[C@H](C/C=C\\CCCCC)O1)O)=C\\C/C=C\\CCCC(=O)O NAME: '
               'hepoxilin B3 REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'C(\\CC)=C\\C/C=C\\C[C@@H](/C=C/C=C\\C/C=C\\CCCC(=O)O)O NAME: '
               '12(S)-HEPE REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/CC(OO)\\C=C\\C=C/CCCC(O)=O NAME: '
               '9-HPETE REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: '
               'CCCCC\\C=C/C[C@@H](O)C(O)\\C=C\\C(O)C\\C=C/CCCC(O)=O NAME: '
               '(5Z,9E,12S,14Z)-8,11,12-trihydroxyicosa-5,9,14-trienoic acid '
               'REASON: WRONGLY CLASSIFIED Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               '[H][C@]12C[C@]3([H])C(=C)C(=O)O[C@@]3([H])[C@@H](O)\\C(C)=C\\CC\\C(C)=C\\CC[C@@]1(CO)O2 '
               'NAME: sinularolide C REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C(O)=O '
               'NAME: (2Z,5Z,8Z,11Z,14Z,17Z)-icosa-2,5,8,11,14,17-hexaenoic '
               'acid REASON: WRONGLY CLASSIFIED Molecule is a nonclassic '
               'icosanoid\n'
               ' * SMILES: O=C(CCC/C=C\\C/C=C\\CC(C(/C=C/C(CCCCC)O)O)O)[O-] '
               'NAME: 11,12,15-trihydroxy-(5Z,8Z,13E)-icosatrienoate REASON: '
               'WRONGLY CLASSIFIED Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'OCCCCC\\C=C/C[C@@H](O)\\C=C\\C=C\\C=C/[C@@H](O)CCCC([O-])=O '
               'NAME: 20-hydroxy-leukotriene B4(1-) REASON: WRONGLY CLASSIFIED '
               'Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(CCC[C@@H](/C=C/C=C\\C/C=C\\C=C\\[C@H](CCCCC)OO)OO)([O-])=O '
               'NAME: (5S,15S)-dihydroperoxy-(6E,8Z,11Z,13E)-icosatetraenoate '
               'REASON: WRONGLY CLASSIFIED Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'O=C1NC[C@H](C=CC=C[C@H]2[C@@H]3[C@@H](C=C1)C=C[C@H]4C[C@H]([C@H]([C@@H]34)O2)O)C '
               'NAME: Cyclamenol B REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC([O-])=O NAME: '
               'arachidonate REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'S([C@@]12N(C(=O)[C@]3(SC)CC=4[C@H](N3C1=O)[C@@H](O)C=CC4)[C@@H]5[C@@H](O)C=CC=C5C2)C '
               'NAME: Boydine A REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: '
               'O=C(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\[C@@H](CCCC)O)[O-] NAME: '
               '16(R)-HETE(1-) REASON: WRONGLY CLASSIFIED Molecule is a '
               'nonclassic icosanoid\n'
               ' * SMILES: O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O '
               'NAME: 5-Oxo-6-trans-leukotriene B4 REASON: WRONGLY CLASSIFIED '
               'Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C1(C)(CO)CCC(C(=C1\\C=C\\C(=C\\C=C\\C(=C\\C(=O)[O-])\\C)\\C)C)=O '
               'NAME: all-trans-4-oxo-16-hydroxyretinoate REASON: WRONGLY '
               'CLASSIFIED Molecule is a nonclassic icosanoid\n'
               ' * SMILES: '
               'C(CCC)C/C=C\\C[C@H]([C@@H](C/C=C\\C/C=C\\CCCC([O-])=O)O)O '
               'NAME: (5Z,8Z,11R,12R,14Z)-11,12-dihydroxyicosatrienoate '
               'REASON: WRONGLY CLASSIFIED Molecule is a nonclassic icosanoid\n'
               ' * SMILES: OC(=O)/C(/CCC=C(C)C)=C\\CC/C(=C/CC/C(/C)=C/CO)/C '
               'NAME: Isonerylgeraniol-18-oic acid REASON: WRONGLY CLASSIFIED '
               'Molecule is a nonclassic icosanoid\n'
               'False negatives: SMILES: '
               'CCCCC[C@H](O)\\C=C\\C=C/C=C/C=C/[C@@H](SC[C@H](N)C(=O)NCC(O)=O)[C@@H](O)CCCC(O)=O '
               'NAME: lipoxin D4 REASON: MISSED Does not have 20 carbons, '
               'which is required for C20 fatty acids\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\C=C/C=C/C=C/[C@@H](SC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(O)=O)[C@@H](O)CCCC(O)=O '
               'NAME: lipoxin C4 REASON: MISSED Does not have 20 carbons, '
               'which is required for C20 fatty acids\n'
               ' * SMILES: '
               'O[C@@H](CCCC(OC)=O)[C@H](O)/C=C/C=C/C=C\\C=C\\[C@@H](O)CCCCC '
               'NAME: 5(S),6(R)-Lipoxin A4 methyl ester REASON: MISSED Does '
               'not have 20 carbons, which is required for C20 fatty acids\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'Insufficient double bonds (at '
                                               'least three required for '
                                               'icosatrienes)'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'Lacking essential functional '
                                               'groups (need at least one OH, '
                                               'epoxy, or COOH group)'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'Does not have at least 20 '
                                               'carbons, which is required for '
                                               'C20 fatty acids'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'Does not have at least 20 '
                                               'carbons, which is required for '
                                               'C20 fatty acids'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'Does not have at least 20 '
                                               'carbons, which is required for '
                                               'C20 fatty acids'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'Insufficient double bonds (at '
                                               'least three required for '
                                               'icosatrienes)'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'Insufficient double bonds (at '
                                               'least three required for '
                                               'icosatrienes)'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'Insufficient double bonds (at '
                                               'least three required for '
                                               'icosatrienes)'},
                                 {   'smiles': 'O=C(N(CCO)CCO)CN1C(=NC=C1)[N+]([O-])=O',
                                     'name': '1H-Imidazole-1-acetamide, '
                                             'N,N-bis(2-hydroxyethyl)-2-nitro-',
                                     'reason': 'Does not have at least 20 '
                                               'carbons, which is required for '
                                               'C20 fatty acids'},
                                 {   'smiles': 'CC1=NN(C(C1)(C2=CN=CC=C2)O)C(=O)CC3=CC=C(C=C3)C',
                                     'name': '1-[5-hydroxy-3-methyl-5-(3-pyridinyl)-4H-pyrazol-1-yl]-2-(4-methylphenyl)ethanone',
                                     'reason': 'Does not have at least 20 '
                                               'carbons, which is required for '
                                               'C20 fatty acids'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'To classify a chemical as a nonclassic icosanoid, I would '
                 'approach this problem by focusing on some key structural '
                 'characteristics:\n'
                 '\n'
                 '1. **C20 Fatty Acid Backbone**: The molecule should consist '
                 'of 20 carbon atoms, typically forming a long carbon chain '
                 'indicative of a fatty acid.\n'
                 '\n'
                 '2. **Double Bonds**: Nonclassic icosanoids usually contain '
                 'multiple double bonds, specifically searching for trends '
                 'like at least three double bonds which characterizes an '
                 'icosatrienoid backbone.\n'
                 '\n'
                 '3. **Functional Groups**: Ensure the presence of essential '
                 'functional groups such as hydroxyl groups (OH), epoxide '
                 'groups, and carboxylic acid groups (COOH) which are '
                 'characteristic of these molecules.\n'
                 '\n'
                 '4. **Exclusion of Classic Derivatives**: The function must '
                 'exclude molecules that match known substructures for classic '
                 'leukotrienes and prostanoids, which are not considered '
                 'nonclassic icosanoids.\n'
                 '\n'
                 'By examining the previous outcomes, the following '
                 'improvements can be made:\n'
                 '\n'
                 '- **Ensure Distinction of C20 Chains**: In cases of false '
                 'negatives, like Lipoxins, the function misunderstood carbon '
                 'counting due to their diversity in branching or cyclic '
                 "forms. It's essential to correctly identify C20 chains even "
                 'if they are embedded in a complex molecular structure.\n'
                 '\n'
                 '- **Enhance Double Bond Detection**: Previous logic filters '
                 'out molecules based on their double bond count, but it needs '
                 'refinement to tolerate different positions and cis-trans '
                 'isomerism in the backbone.\n'
                 '\n'
                 '- **Functional Group Identification**: Better identification '
                 'logic for crucial functional groups, ensuring all spatial '
                 'isomers are accounted.\n'
                 '\n'
                 '- **Classic icosanoid Exclusion**: Ensure excluded patterns '
                 'are specific and distinguish correctly versus the structural '
                 'diversity in nonclassic icosanoids.\n'
                 '\n'
                 "Now, let's write the improved definition of "
                 '`is_nonclassic_icosanoid`:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 230,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.7198879551820728,
    'negative_predictive_value': 1.0}