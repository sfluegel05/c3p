"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring nitrogen-containing compound (mostly heterocyclic),
    occurring mostly in plants, but also found in bacteria, fungi, and animals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    # Check for nitrogen atoms in rings (heterocyclic nitrogen)
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    has_heterocyclic_nitrogen = False
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_atomic_nums = set(atom.GetAtomicNum() for atom in ring_atoms)
        if 7 in ring_atomic_nums:
            has_heterocyclic_nitrogen = True
            break

    if not has_heterocyclic_nitrogen:
        return False, "No nitrogen atoms found in rings (no heterocyclic nitrogen)"

    # Exclude amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics
    # Check for standard amino acid backbone
    amino_acid = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")
    if mol.HasSubstructMatch(amino_acid):
        return False, "Structure matches amino acid backbone"

    # Check for peptide bond pattern (N-C(=O)-C)
    peptide_bond = Chem.MolFromSmarts("N-C(=O)-C")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond)
    if len(peptide_bonds) > 3:  # Arbitrary threshold to exclude peptides/proteins
        return False, "Contains multiple peptide bonds (possible peptide or protein)"

    # Check for nucleic acid bases
    purine = Chem.MolFromSmarts("c1ncnc2ncnc12")
    pyrimidine = Chem.MolFromSmarts("c1ccncn1")
    if mol.HasSubstructMatch(purine) or mol.HasSubstructMatch(pyrimidine):
        return False, "Contains purine or pyrimidine base (possible nucleotide or nucleic acid)"

    # Check for sugar moieties (e.g., furanose or pyranose rings)
    sugar_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1O")
    if mol.HasSubstructMatch(sugar_pattern):
        return False, "Contains sugar moiety (possible amino sugar or glycoside)"

    # Check for antibiotics common substructures (e.g., beta-lactam)
    beta_lactam = Chem.MolFromSmarts("C1C(=O)NC1")
    if mol.HasSubstructMatch(beta_lactam):
        return False, "Contains beta-lactam ring (possible antibiotic)"

    # If all checks pass, classify as alkaloid
    return True, "Molecule contains heterocyclic nitrogen characteristic of alkaloids"

__metadata__ = {   'chemical_class': {   'name': 'alkaloid',
                              'definition': 'Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals. By extension, certain neutral compounds biogenetically related to basic alkaloids are also classed as alkaloids. Amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars and antibiotics are not normally regarded as alkaloids. Compounds in which the nitrogen is  exocyclic (dopamine, mescaline, serotonin, etc.) are usually classed as amines rather than alkaloids.',
                              'parents': []},
    'message': None,
    'success': True}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22315',
                          'name': 'alkaloid',
                          'definition': 'Any of the naturally occurring, basic '
                                        'nitrogen compounds (mostly '
                                        'heterocyclic) occurring mostly in the '
                                        'plant kingdom, but also found in '
                                        'bacteria, fungi, and animals. By '
                                        'extension, certain neutral compounds '
                                        'biogenetically related to basic '
                                        'alkaloids are also classed as '
                                        'alkaloids. Amino acids, peptides, '
                                        'proteins, nucleotides, nucleic acids, '
                                        'amino sugars and antibiotics are not '
                                        'normally regarded as alkaloids. '
                                        'Compounds in which the nitrogen is  '
                                        'exocyclic (dopamine, mescaline, '
                                        'serotonin, etc.) are usually classed '
                                        'as amines rather than alkaloids.',
                          'parents': ['CHEBI:35352'],
                          'xrefs': ['Wikipedia:Alkaloid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6710526315789473 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[H]\\C(C)=C1/CN(C)[C@]2([H])Cc3c(CC[C@]1([H])[C@]2([H])C)[nH]c1ccccc31 '
               'NAME: vobasan REASON: CORRECT Molecule contains basic nitrogen '
               'in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'C1OCC=2[C@@]3([C@H](OC/C(/C[C@H]([C@@H]1C)C)=C(/[H])\\C)CCN3CC2)[H] '
               'NAME: senecionan REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'O[C@H]1[C@@]([C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(C1)[H])C)([C@@H]([C@]5(NC[C@H](CC5)C)[H])C)[H] '
               'NAME: Teinemine REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'O([C@@H]1C2N(C(C1)C[C@H](OC(=O)/C(/C)=C/C)C2)C)C(=O)/C(/C)=C/C '
               'NAME: L-3alpha,6beta-Ditigloyloxytropane REASON: CORRECT '
               'Molecule contains basic nitrogen in a heterocyclic ring '
               'characteristic of alkaloids\n'
               ' * SMILES: COc1cc2CCN(C)[C@@]3(Cc4ccc5OCOc5c4C3=O)c2cc1O NAME: '
               'parfumine REASON: CORRECT Molecule contains basic nitrogen in '
               'a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'O([C@@H](C[C@]12[C@]3([C@@]4(N5[C@]1([C@H](CC[C@]2([C@@](C5)(CC4)[H])C)C(C)C)[H])CCC3)[H])C(=O)[C@]6([C@]7(O[C@@](OC6)(CC7)C)[H])C)C(=O)C '
               'NAME: Daphniphylline REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               '[H][C@@]1(CN2CC[C@@]34C(Nc5ccccc35)=C(C(=O)OC)[C@@]1([H])C[C@]24[H])C(C)=O '
               'NAME: (-)-Alstolucine F REASON: CORRECT Molecule contains '
               'basic nitrogen in a heterocyclic ring characteristic of '
               'alkaloids\n'
               ' * SMILES: CC[C@@H](C)C[C@@H]1CCCN(CC)C1 NAME: '
               "(2'R,3S)-stenusine REASON: CORRECT Molecule contains basic "
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'COC(=O)C1=CO[C@@H](C)[C@H]2CN3CCc4c([nH]c5ccc(OC)cc45)[C@@H]3C[C@H]12 '
               'NAME: Cabucine REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               '[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@@]6(O[C@]5(C1)[H])NC[C@H](C)CC6)C)[H])[H])[H])[H] '
               'NAME: solasodine REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: C1CCN2C[C@@H]3C[C@H](C=C4CCCN[C@@H]34)[C@H]2C1 '
               'NAME: Aloperine REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CN1C2=C(C=CC(=C2)OC)C3=C1[C@@H](N(CC34CCN(CC4)S(=O)(=O)C)CC5=CC(=CC=C5)OC)CO '
               'NAME: '
               "[(1R)-7-methoxy-2-[(3-methoxyphenyl)methyl]-9-methyl-1'-methylsulfonyl-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]methanol "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CCS(=O)(=O)N1CCC2(CC1)CN([C@@H](C3=C2C4=C(N3)C=C(C=C4)OC)CO)S(=O)(=O)C5=CC=C(C=C5)C '
               'NAME: '
               "[(1S)-1'-ethylsulfonyl-7-methoxy-2-(4-methylphenyl)sulfonyl-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1 '
               'NAME: tropan-3alpha-yl 3-hydroxy-2-phenylpropanoate REASON: '
               'CORRECT Molecule contains basic nitrogen in a heterocyclic '
               'ring characteristic of alkaloids\n'
               ' * SMILES: O=C(O)[C@H]1N[C@@H](C=2NC3=C(C2C1)C=CC(=C3)O)C '
               'NAME: Brunnein B REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CC(C)NC(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)CC5=CC=C(C=C5)Cl '
               'NAME: '
               "(1R)-2-[(4-chlorophenyl)methyl]-1-(hydroxymethyl)-7-methoxy-9-methyl-N-propan-2-yl-1'-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]carboxamide "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](NCC34CCN(CC4)CC5CCCC5)CO NAME: '
               "[(1S)-1'-(cyclopentylmethyl)-7-methoxy-9-methyl-1-spiro[2,3-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CC1=CC=C(C=C1)S(=O)(=O)N2CC3(CCN(CC3)CC4=NC=CS4)C5=C([C@@H]2CO)NC6=C5C=CC(=C6)OC '
               'NAME: '
               "[(1R)-7-methoxy-2-(4-methylphenyl)sulfonyl-1'-(2-thiazolylmethyl)-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'C\\C=C1\\CN2[C@H]3Cc4c([C@@H]2C[C@@H]1C3CO)n(C)c1ccccc41 NAME: '
               'Affinisine REASON: CORRECT Molecule contains basic nitrogen in '
               'a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'C12=C(C=CC=C1)NC=3[C@]([C@]4(/C(/CN(CC4)CCC32)=C\\C)[H])(C(=O)OC)CO '
               'NAME: 15beta-stemmadenine REASON: CORRECT Molecule contains '
               'basic nitrogen in a heterocyclic ring characteristic of '
               'alkaloids\n'
               ' * SMILES: '
               'O1[C@@]2([C@]34[C@](O)([C@](N(CC3)C)(CC5=C4C1=C(O)C=C5)[H])C=C([C@H]6N([C@]7(C=8[C@]9(C6)[C@@](OC=%10C9=C(C7)C=CC%10OC)(C(OC)=CC8)[H])[H])C)C2=O)[H] '
               'NAME: Somniferine REASON: CORRECT Molecule contains basic '
               'nitrogen in a heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'COC1=CC2=C(C=C1)C3=C(N2)[C@@H](N(CC34CN(C4)CC5CCOCC5)C(=O)NC6=CC(=CC=C6)F)CO '
               'NAME: '
               "(1R)-N-(3-fluorophenyl)-1-(hydroxymethyl)-7-methoxy-1'-(4-oxanylmethyl)-2-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]carboxamide "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CCS(=O)(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)C(=O)C5=CC=C(C=C5)F '
               'NAME: '
               "[(1R)-1'-ethylsulfonyl-1-(hydroxymethyl)-7-methoxy-9-methyl-2-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]-(4-fluorophenyl)methanone "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'CC[C@H]1CN2CCC3=CC(=C(C=C3C2C[C@@H]1CC4C5=CC(=C(C=C5CCN4)OC)OC)OC)OC '
               'NAME: '
               '(2S,3R)-2-[(6,7-dimethoxy-1,2,3,4-tetrahydroisoquinolin-1-yl)methyl]-3-ethyl-9,10-dimethoxy-2,3,4,6,7,11b-hexahydro-1H-benzo[a]quinolizine '
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               ' * SMILES: '
               'COC1=CC2=C(C=C1)C3=C(N2)[C@@H](N(CC34CN(C4)C(=O)C5=CC=C(C=C5)F)CC6CC6)CO '
               'NAME: '
               "[(1R)-2-(cyclopropylmethyl)-1-(hydroxymethyl)-7-methoxy-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]yl]-(4-fluorophenyl)methanone "
               'REASON: CORRECT Molecule contains basic nitrogen in a '
               'heterocyclic ring characteristic of alkaloids\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CCC(=O)N1CC2(CCN(CC2)C(=O)CC3=CC=CC=C3)C4=C([C@H]1CO)N(C5=C4C=CC(=C5)OC)C '
               'NAME: '
               "1-[(1S)-1-(hydroxymethyl)-7-methoxy-9-methyl-1'-(1-oxo-2-phenylethyl)-2-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]-1-propanone "
               'REASON: MISSED No basic (sp3-hybridized) nitrogen atoms found\n'
               ' * SMILES: '
               'CC1=CC=C(C=C1)S(=O)(=O)N2CC3(CCN(CC3)C(=O)CC4=CN=CC=C4)C5=C([C@H]2CO)NC6=C5C=CC(=C6)OC '
               'NAME: '
               "1-[(1S)-1-(hydroxymethyl)-7-methoxy-2-(4-methylphenyl)sulfonyl-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]-2-(3-pyridinyl)ethanone "
               'REASON: MISSED Contains peptide bonds (possible peptide or '
               'protein)\n'
               ' * SMILES: '
               'O=C1N2[C@@]3([C@@]4([C@]5(N(CC4)C/C(/[C@](C5)(C3=CC1)[H])=C/CO)[H])C=6C2=CC(=C(C6)OC)OC)[H] '
               'NAME: isobrucine REASON: MISSED Contains peptide bonds '
               '(possible peptide or protein)\n'
               ' * SMILES: '
               'O[C@H]1[C@]2(N(C[C@@]3([C@@]4(N(CCC[C@@]42[H])CCC3)[H])[H])C(=O)C=C1)[H] '
               'NAME: (-)-12beta-Hydroxylsophocarpine REASON: MISSED Contains '
               'peptide bonds (possible peptide or protein)\n'
               ' * SMILES: [O-][N+](=O)c1c(Cl)cccc1-c1c[nH]cc1Cl NAME: '
               'pyrrolnitrin REASON: MISSED No basic (sp3-hybridized) nitrogen '
               'atoms found\n'
               ' * SMILES: COC(=O)\\C=C/NC(=O)c1cc2c3ccccc3[nH]c2c(n1)C(C)=O '
               'NAME: Stellarine B REASON: MISSED No basic (sp3-hybridized) '
               'nitrogen atoms found\n'
               ' * SMILES: '
               'CN1C2=C(C=CC(=C2)OC)C3=C1[C@@H](N(CC34CN(C4)C(=O)CC5CC5)CC6=CC=CC=C6Cl)CO '
               'NAME: '
               "1-[(1R)-2-[(2-chlorophenyl)methyl]-1-(hydroxymethyl)-7-methoxy-9-methyl-1'-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]-2-cyclopropylethanone "
               'REASON: MISSED Contains peptide bonds (possible peptide or '
               'protein)\n'
               ' * SMILES: '
               '[H][C@]12[C@@H](O)C[C@H](N)C(=O)N1CC1=C(C=CC(O)=C1O)[C@H]2O '
               'NAME: fumisoquin A REASON: MISSED Contains peptide bonds '
               '(possible peptide or protein)\n'
               ' * SMILES: '
               'CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](NCC34CN(C4)C(=O)CC5CC5)CO NAME: '
               "2-cyclopropyl-1-[(1S)-1-(hydroxymethyl)-7-methoxy-9-methyl-1'-spiro[2,3-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]yl]ethanone "
               'REASON: MISSED Contains peptide bonds (possible peptide or '
               'protein)\n'
               ' * SMILES: '
               '[C@@H]1([C@H]([C@@]2([C@@H](OC(=O)C)[C@H]([C@]3([C@](O)([C@]42[C@H](OC(=O)C)[C@]1([H])[C@@](COC(C5=C(CC[C@@H](C(O3)=O)C)N=CC=C5)=O)(C)O4)C)[H])O)COC(=O)C)OC(=O)C)OC(C)=O '
               'NAME: 2-O-deacetyleuonine REASON: MISSED No basic '
               '(sp3-hybridized) nitrogen atoms found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'C[C@H](Cc1ccccc1)N(C)CC#C',
                                      'name': '(-)-selegiline',
                                      'reason': 'No nitrogen atoms found in '
                                                'rings (no heterocyclic '
                                                'nitrogen)'},
                                  {   'smiles': 'OC(=O)[C@@H]1CC(=C/C=O)\\C=C(N1)C(O)=O',
                                      'name': 'betalamic acid',
                                      'reason': 'Structure matches amino acid '
                                                'backbone'},
                                  {   'smiles': 'O=C(O)[C@H]1N[C@@H](C=2NC3=C(C2C1)C=CC(=C3)O)C',
                                      'name': 'Brunnein B',
                                      'reason': 'Structure matches amino acid '
                                                'backbone'},
                                  {   'smiles': 'CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](N(CC34CCN(CC4)CC5=CN=CN=C5)S(=O)(=O)C)CO',
                                      'name': "[(1S)-7-methoxy-9-methyl-2-methylsulfonyl-1'-(5-pyrimidinylmethyl)-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                      'reason': 'Contains purine or pyrimidine '
                                                'base (possible nucleotide or '
                                                'nucleic acid)'},
                                  {   'smiles': 'O=C1C(=C(C=2C3=C(C=CC=C3)NC2)C(C14N[C@H](C(=O)O)CC5=C4NC=6C=CC=CC56)=O)C=7C8=C(C=CC=C8)NC7',
                                      'name': 'Pityriarubin A',
                                      'reason': 'Structure matches amino acid '
                                                'backbone'},
                                  {   'smiles': 'COc1ccc2C(=O)c3c(OC)cc(CCN(C)C)c(O)c3-c2c1O',
                                      'name': 'caulophine',
                                      'reason': 'No nitrogen atoms found in '
                                                'rings (no heterocyclic '
                                                'nitrogen)'},
                                  {   'smiles': 'OC(=O)C1NC(C=2NC=3C(C2C1)=CC=CC3)C',
                                      'name': '(1xi,3xi)-1,2,3,4-Tetrahydro-1-methyl-beta-carboline-3-carboxylic '
                                              'acid',
                                      'reason': 'Structure matches amino acid '
                                                'backbone'},
                                  {   'smiles': '[H][C@]12N[C@@H](C)C(=O)N1c1ccccc1[C@@]21C[C@@H]2C(=O)N[C@@](C)(O1)c1nc3ccccc3c(=O)n21',
                                      'name': 'fumiquinazoline C',
                                      'reason': 'Contains purine or pyrimidine '
                                                'base (possible nucleotide or '
                                                'nucleic acid)'},
                                  {   'smiles': 'CNC[C@H](O)c1ccc(O)cc1',
                                      'name': 'D-synephrine',
                                      'reason': 'No nitrogen atoms found in '
                                                'rings (no heterocyclic '
                                                'nitrogen)'},
                                  {   'smiles': 'C/C=C\\1/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@]36OC(C7=C4C=CC(=C7O)O)=O)[H])[H])(C(=O)OC)[H]',
                                      'name': 'Voacalgine A',
                                      'reason': 'Structure matches amino acid '
                                                'backbone'}],
    'attempt': 1,
    'reasoning': 'Upon analyzing the previous program and the outcomes, it '
                 'seems that the code was too restrictive in its criteria for '
                 'classifying alkaloids. Specifically:\n'
                 '\n'
                 '1. **Requirement for basic (sp3-hybridized) nitrogen '
                 'atoms**: The previous code only considered molecules with '
                 'basic, sp3-hybridized nitrogen atoms as alkaloids. However, '
                 'the definition includes certain neutral compounds '
                 'biogenetically related to basic alkaloids. Therefore, '
                 'excluding molecules with non-basic nitrogen atoms might lead '
                 'to false negatives.\n'
                 '\n'
                 '2. **Exclusion of molecules containing peptide bonds**: The '
                 'code excluded any molecule containing peptide bonds, but '
                 'some alkaloids might have amide linkages or peptide-like '
                 'structures. This broad exclusion can wrongly filter out true '
                 'alkaloids.\n'
                 '\n'
                 '3. **Exclusion based on exocyclic nitrogen**: The code '
                 'excluded compounds with exocyclic nitrogen atoms if they '
                 "didn't also have nitrogen in a heterocyclic ring. Some "
                 'alkaloids may have exocyclic nitrogen without fitting this '
                 'strict criterion.\n'
                 '\n'
                 'To improve the classification:\n'
                 '\n'
                 '- **Relax the requirement on nitrogen hybridization**: '
                 'Consider any nitrogen atom in a ring (heterocyclic '
                 'nitrogen), regardless of its hybridization or basicity.\n'
                 '\n'
                 '- **Allow for neutral nitrogen compounds**: Accept compounds '
                 'that are biogenetically related to alkaloids, even if they '
                 "don't have basic nitrogen atoms.\n"
                 '\n'
                 '- **Refine exclusion criteria**: Instead of broadly '
                 'excluding all molecules with peptide bonds, focus on '
                 'excluding standard amino acids, peptides, proteins, '
                 'nucleotides, nucleic acids, amino sugars, and antibiotics '
                 'more specifically.\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1563,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9398677089597114,
    'f1': 0.9690018598884067,
    'accuracy': 0.9398677089597114,
    'negative_predictive_value': 0.0}