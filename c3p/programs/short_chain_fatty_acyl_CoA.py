"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    Short-chain fatty acids typically have 2-6 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA core structure
    coa_pattern = Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]2:[#6](:[#7]:[#6]:[#7]:[#6]:2:1)[#7]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA adenine base structure"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester linkage"

    # Check for pantetheine arm
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm"

    # Find the acyl portion by first identifying the thioester carbon
    thioester_carbon = thioester_matches[0][0]
    
    # Create visited atoms set starting from thioester carbon
    visited = set()
    queue = [thioester_carbon]
    acyl_carbons = 0
    
    # Breadth-first search to count carbons in acyl chain
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        
        atom = mol.GetAtomWithIdx(current)
        # Count carbon atoms (excluding the thioester carbon)
        if atom.GetAtomicNum() == 6 and current != thioester_carbon:
            acyl_carbons += 1
            
        # Add neighbors to queue, but stop at sulfur (thioester linkage)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 16:  # not sulfur
                queue.append(neighbor.GetIdx())

    if acyl_carbons < 2:
        return False, f"Acyl chain too short ({acyl_carbons} carbons)"
    if acyl_carbons > 6:
        return False, f"Acyl chain too long ({acyl_carbons} carbons) for short-chain fatty acid"

    # Check for common modifications in the acyl portion
    modifications = []
    mod_patterns = {
        "hydroxy": ("[OX2H1]", "hydroxyl"),
        "amino": ("[NX3H2]", "amino"),
        "unsaturated": ("[CX3]=[CX3]", "unsaturated"),
        "methyl-branched": ("[CH3][CH]([!H])[!H]", "methyl-branched")
    }
    
    for mod_name, (pattern, desc) in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            modifications.append(desc)

    mod_text = " with " + ", ".join(modifications) if modifications else ""
    return True, f"Short-chain ({acyl_carbons} carbons) fatty acyl-CoA{mod_text}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61905',
                          'name': 'short-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any short-chain fatty acid.',
                          'parents': ['CHEBI:37554'],
                          'xrefs': [   'PMID:20371198',
                                       'PMID:20376488',
                                       'PMID:20429031'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC([C@H](CC)C)=O)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (S)-2-methylbutanoyl-CoA REASON: MISSED Acyl chain too '
               'long (26 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)[C@@H](CC)O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (R)-2-hydroxybutanoyl-CoA REASON: MISSED Acyl chain too '
               'long (25 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'C[C@H](N)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: L-3-aminobutanoyl-CoA REASON: MISSED Acyl chain too long '
               '(25 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'C[C@H](O)[C@H](C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2S,3S)-3-hydroxy-2-methylbutanoyl-CoA REASON: MISSED '
               'Acyl chain too long (26 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'CC(=C)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-methylbut-3-enoyl-CoA REASON: MISSED Acyl chain too '
               'long (26 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (R)-3-hydroxybutanoyl-CoA REASON: MISSED Acyl chain too '
               'long (25 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'CC(C)C=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 4-methylpent-2-enoyl-CoA REASON: MISSED Acyl chain too '
               'long (27 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: butyryl-CoA REASON: MISSED Acyl chain too long (25 '
               'carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C(CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: 2-hydroxybutanoyl-CoA REASON: MISSED Acyl chain too long '
               '(25 carbons) for short-chain fatty acid\n'
               ' * SMILES: '
               'CC(C)(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxyisovaleryl-CoA REASON: MISSED Acyl chain too '
               'long (26 carbons) for short-chain fatty acid\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Missing thioester linkage'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Missing CoA adenine base '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'COc1cc(ccc1O)C(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12',
                                      'name': '3-hydroxy-3-(4-hydroxy-3-methoxyphenyl)propanoyl-CoA',
                                      'reason': 'Acyl chain too long (9 '
                                                'carbons) for short-chain '
                                                'fatty acid'}],
    'attempt': 3,
    'reasoning': 'Let me analyze what went wrong in the previous attempt and '
                 'suggest improvements:\n'
                 '\n'
                 '1. Main Issue: Carbon Counting\n'
                 'The program was incorrectly counting carbons in the acyl '
                 'chain. When it tried to isolate the acyl fragment, it was '
                 'likely including parts of the CoA structure, leading to '
                 'incorrect carbon counts and false negatives.\n'
                 '\n'
                 '2. Fragment Analysis Problems\n'
                 'The fragmentation approach used might not have cleanly '
                 'separated the acyl chain from the CoA portion, causing '
                 'miscalculation of chain length.\n'
                 '\n'
                 '3. Missing Validation\n'
                 'Some compounds were rejected because the program thought '
                 'they had too many carbons, but it was actually counting the '
                 'entire molecule instead of just the fatty acid portion.\n'
                 '\n'
                 "Here's an improved version that should fix these issues:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 29,
    'num_false_positives': 100,
    'num_true_negatives': 67449,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.2248062015503876,
    'recall': 0.9666666666666667,
    'f1': 0.36477987421383645,
    'accuracy': 0.9985054528773731,
    'negative_predictive_value': 0.9999851742031134}