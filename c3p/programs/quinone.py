"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:26421 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclic ketone patterns
    # This pattern matches both para and ortho quinones
    para_quinone = Chem.MolFromSmarts("[#6;R]=O.[#6;R]=O")
    ortho_quinone = Chem.MolFromSmarts("[#6;R](=O)-[#6;R]-[#6;R](=O)")
    
    has_para = mol.HasSubstructMatch(para_quinone)
    has_ortho = mol.HasSubstructMatch(ortho_quinone)
    
    if not (has_para or has_ortho):
        return False, "Must contain cyclic ketone groups in quinone arrangement"

    # Count ring systems containing ketones
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(ketone_matches) < 2:
        return False, "Must have at least two ketone groups"
    
    # Check if ketones are in conjugated system
    ring_info = mol.GetRingInfo()
    ring_systems = ring_info.AtomRings()
    
    # Get the atoms between ketones
    connecting_atoms = []
    for ring in ring_systems:
        ketones_in_ring = []
        for ketone in ketone_matches:
            if ketone[0] in ring:
                ketones_in_ring.append(ketone[0])
        if len(ketones_in_ring) >= 2:
            connecting_atoms.extend(ring)
    
    if not connecting_atoms:
        return False, "Ketone groups must be in same ring system"
    
    # Check for conjugation or aromaticity in connecting atoms
    has_conjugation = False
    for atom_idx in connecting_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            has_conjugation = True
            break
        for bond in atom.GetBonds():
            if bond.GetIsConjugated():
                has_conjugation = True
                break
    
    if not has_conjugation:
        return False, "Must have conjugated system between ketones"
    
    # Additional check for polycyclic systems
    polycyclic = len(ring_systems) > 1
    
    # Construct reason string
    if has_para:
        structure_type = "para"
    else:
        structure_type = "ortho"
    
    system_type = "polycyclic" if polycyclic else "monocyclic"
    reason = f"Contains {structure_type}-quinone structure in {system_type} conjugated system"
    
    return True, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36141',
                          'name': 'quinone',
                          'definition': 'Compounds having a fully conjugated '
                                        'cyclic dione structure, such as that '
                                        'of benzoquinones, derived from '
                                        'aromatic compounds by conversion of '
                                        'an even number of -CH= groups into '
                                        '-C(=O)- groups with any necessary '
                                        'rearrangement of double bonds '
                                        '(polycyclic and heterocyclic '
                                        'analogues are included).',
                          'parents': ['CHEBI:3992'],
                          'xrefs': ['Wikipedia:Quinone'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: O=C1C=CC(=O)c2cc3cc4ccccc4cc3cc12 '
               'NAME: tetracene-1,4-dione REASON: MISSED Must have '
               'characteristic quinone pattern\n'
               ' * SMILES: ClC1=C2C(=O)C3=C(C(OC)=CC(=C3)C)C(C2=C(O)C=C1O)=O '
               'NAME: 5-Chloro-6,8-dihydroxy-1-methoxy-3-methyanthraquinone '
               'REASON: MISSED Must have at least 2 cyclic ketone groups\n'
               ' * SMILES: CC(O)c1cc2C(=O)c3cccc(O)c3C(=O)c2o1 NAME: '
               "8-hydroxy-2-(1'-hydroxyethyl)naphtho-[2,3-b]-furan-4,9-dione "
               'REASON: MISSED Must have at least 2 cyclic ketone groups\n'
               ' * SMILES: '
               'O=C1C2=C(OC)C=C(COC(=O)C)C=C2C(=O)C=3C1=C(O)C=C(O)C3 NAME: '
               'omega-acetylcarviolin REASON: MISSED Must have at least 2 '
               'cyclic ketone groups\n'
               ' * SMILES: COC1=C(OC)C(=O)C(CCCCCCCCCCO)=C(C)C1=O NAME: '
               'idebenone REASON: MISSED Must have characteristic quinone '
               'pattern\n'
               ' * SMILES: '
               'O=C1C2=C(O)C=C(C)C=C2C(=O)C3=C1C=4O[C@]56O[C@@H](C)[C@H]7C6(C4C=C3)C(=C(O)C8=C5C=CC=C8O)C(=O)C=C7C '
               'NAME: '
               '(1R,18R,31S)-7,23,25-trihydroxy-9,19,31-trimethyl-2,30-dioxaoctacyclo[16.11.2.01,17.03,16.04,13.06,11.017,22.024,29]hentriaconta-3(16),4(13),6(11),7,9,14,19,22,24(29),25,27-undecaene-5,12,21-trione '
               'REASON: MISSED Ketone groups must be conjugated\n'
               ' * SMILES: '
               'O=C1C2=C(O)C=3C(=C4C(O)=CC(CCCCC)=CC4=CC3)C=C2C(=O)C=5C1=C(O)C=C(O)C5 '
               'NAME: Bequinostatin D REASON: MISSED Must have at least 2 '
               'cyclic ketone groups\n'
               ' * SMILES: COc1c(O)c(C)cc2C(=O)c3cccc(O)c3C(=O)c12 NAME: '
               'Obtusifolin REASON: MISSED Must have at least 2 cyclic ketone '
               'groups\n'
               ' * SMILES: '
               'C1=CC=CC2=C1C(C3=C(C2=O)C(=C(C(=C3N)O)S(O)(=O)=O)O)=O NAME: '
               'nuclear fast red free acid REASON: MISSED Must have at least 2 '
               'cyclic ketone groups\n'
               ' * SMILES: '
               'O=C1C2=C(O)C(=CC=C2C(=O)C3=C1C(=CC(=C3)O)C)[C@@H]4O[C@@H]([C@@H](O)[C@@H](C4)N(C)C)C '
               'NAME: Lactoquinomycin E REASON: MISSED Must have at least 2 '
               'cyclic ketone groups\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1=CC=CC(=C1)CC[C@H](N[C@H](C(N2C(N(C[C@H]2C(O)=O)C)=O)=O)C)C(=O)O',
                                     'name': 'imidaprilat',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)CCCN2C=C(CO[C@H]1CN(C)CC3=CC(=CC=C3)O)N=N2)[C@@H](C)CO',
                                     'name': '(8S,9R)-9-[[(3-hydroxyphenyl)methyl-methylamino]methyl]-6-[(2S)-1-hydroxypropan-2-yl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'ClC1=C2C(C3=C(C(=O)CC)C(O[C@@]3(C1=O)C)=O)=CN(CC=C(C)C)C(=C2)/C=C/[C@H](CC)C',
                                     'name': '4′,5′-dinor-5′-deoxy-N-(3-methyl-2-butenyl)-2-aza-2-deoxychaetoviridin '
                                             'A',
                                     'reason': 'Ketone groups must be in same '
                                               'ring system'},
                                 {   'smiles': 'O1[C@@H](O[C@H]2[C@@H](O)[C@H](O[C@H](O)[C@@H]2O)CO)[C@H](NC(=O)C)[C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H](O)[C@H]1CO',
                                     'name': 'Lacto-N-triose I',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'O=C1OCC(=C1/C=C\\C)CC[C@]2(OC(=O)[C@@H](C2)O)C',
                                     'name': 'Pezizolide D',
                                     'reason': 'Ketone groups must be in same '
                                               'ring system'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)C2=C(N=CC(=C2)C3=CCCCC3)O[C@H]1CN(C)S(=O)(=O)C4=CN(C=N4)C)[C@H](C)CO',
                                     'name': 'N-[[(2R,3S)-8-(1-cyclohexenyl)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-2-yl]methyl]-N,1-dimethyl-4-imidazolesulfonamide',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CC)C(O)=O)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@@H](CO)O[C@@]1(C[C@H](O)[C@@H](NC(=O)CC)[C@@]([H])(O1)[C@H](O)[C@H](O)CO)C(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O)C(O)=O',
                                     'name': '[8)-alpha-Neu5Pr-(2->]10',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'O=C1C2=C(OC[C@@]1(O)C3=CC=C(O)C=C3)C=C(O)C=C2O',
                                     'name': 'Dothideoisoflavanone',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'C1CN(CCN1C2=C3C=CC=NC3=C(C=C2)[N+](=O)[O-])C(=O)C4=CC=C(C=C4)Cl',
                                     'name': '(4-chlorophenyl)-[4-(8-nitro-5-quinolinyl)-1-piperazinyl]methanone',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone arrangement'},
                                 {   'smiles': 'O=C(N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)[C@@H](NC(=O)CN)CCC(O)=O',
                                     'name': 'Gly-Glu-Tyr',
                                     'reason': 'Must contain cyclic ketone '
                                               'groups in quinone '
                                               'arrangement'}],
    'sample_false_negatives': [   {   'smiles': 'ClC1=C(OC)C(C2=C3C(O)=C4C(=O)C5=C(O)C=C(O)C=C5C(C4(CC(=O)C)C(C3=CC(=C2)OC)=O)(C)C)=C(C)C=C1OC',
                                      'name': 'Naphthacemycin C2',
                                      'reason': 'Ketone groups must be in same '
                                                'ring system'},
                                  {   'smiles': 'ClC1=C2C(C(O)=C3C(=O)C=4C(O)=CC(=CC4C(C3(C2=O)CC(=O)C)(C)C)O)=C(C5=C(OC)C=C(O)C=C5C)C=C1OC',
                                      'name': 'Naphthacemycin C1',
                                      'reason': 'Ketone groups must be in same '
                                                'ring system'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous code and the failures:\n'
                 '\n'
                 '1. The main issue seems to be that the code is too '
                 'restrictive, particularly in the quinone pattern matching. '
                 'Many valid quinones are failing the checks.\n'
                 '\n'
                 '2. Key observations from failed examples:\n'
                 '- Many valid quinones like tetracene-1,4-dione and idebenone '
                 'were rejected due to the overly specific quinone pattern\n'
                 '- The conjugation check might be too strict\n'
                 "- The code doesn't properly handle different types of "
                 'quinone structures (p-quinones, o-quinones, etc.)\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '- Simplify the core quinone detection to focus on the '
                 'essential features: cyclic structure with conjugated '
                 'ketones\n'
                 '- Better handling of aromatic systems\n'
                 '- More flexible pattern matching for different quinone '
                 'types\n'
                 '- Remove overly restrictive patterns\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1155,
    'num_false_positives': 100,
    'num_true_negatives': 1495,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9203187250996016,
    'recall': 0.9982713915298185,
    'f1': 0.9577114427860697,
    'accuracy': 0.9629360465116279,
    'negative_predictive_value': 0.9986639946559787}