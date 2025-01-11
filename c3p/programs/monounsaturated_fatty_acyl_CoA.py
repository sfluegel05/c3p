"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_monounsaturated_fatty_acyl_CoA, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety components
    # Look for adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for phosphate groups (need at least 3 for CoA)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups for CoA (found {phosphate_matches}, need ≥3)"

    # Look for thioester linkage (more flexible pattern)
    thioester_pattern = Chem.MolFromSmarts("[S;X2][C;X3](=[O;X1])")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Look for amide bonds characteristic of pantetheine
    amide_pattern = Chem.MolFromSmarts("[NX3][C;X3](=[O;X1])[CX4]")
    amide_matches = len(mol.GetSubstructMatches(amide_pattern))
    if amide_matches < 2:
        return False, f"Missing pantetheine amide bonds (found {amide_matches}, need ≥2)"

    # Count non-aromatic double bonds
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Filter out aromatic bonds
    non_aromatic_double_bonds = []
    for bond in double_bonds:
        atom1, atom2 = mol.GetAtomWithIdx(bond[0]), mol.GetAtomWithIdx(bond[1])
        if not (atom1.IsInRing() and atom2.IsInRing() and atom1.GetIsAromatic() and atom2.GetIsAromatic()):
            non_aromatic_double_bonds.append(bond)
    
    if len(non_aromatic_double_bonds) != 1:
        return False, f"Found {len(non_aromatic_double_bonds)} C=C double bonds, need exactly 1"

    # Verify the presence of a fatty acid chain
    # Get carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:  # CoA itself has about 23 carbons, so fatty acyl-CoA should have more
        return False, f"Carbon count ({carbon_count}) too low for fatty acyl-CoA"

    # Verify the double bond is in the fatty acyl portion
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        sulfur_idx = thioester_matches[0][0]
        for db_atoms in non_aromatic_double_bonds:
            # Check if we can find a path from the sulfur to both double bond carbons
            # that doesn't go through the CoA part
            paths_exist = any(len(Chem.FindAllPathsOfLengthN(mol, 2, sulfur_idx, db_atom)) > 0 
                            for db_atom in db_atoms)
            if not paths_exist:
                return False, "Double bond not in fatty acyl chain"

    return True, "Contains CoA moiety and single double bond in fatty acyl chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139575',
                          'name': 'monounsaturated fatty acyl-CoA',
                          'definition': 'Any unsaturated fatty acyl-CoA in '
                                        'which the fatty acyl chain contains '
                                        'one carbon-carbon double bond.',
                          'parents': ['CHEBI:51006'],
                          'xrefs': ['PMID:17277381', 'PMID:3382682'],
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
               'CCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-tetradec-2-enoyl-CoA REASON: MISSED No CoA moiety '
               'found (missing pantetheine)\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z)-octadecenoyl-CoA REASON: MISSED No CoA moiety '
               'found (missing pantetheine)\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC=C '
               'NAME: vinylacetyl-CoA REASON: MISSED No CoA moiety found '
               '(missing pantetheine)\n'
               ' * SMILES: '
               'C\\C=C\\CCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (14E)-hexadecenoyl-CoA REASON: MISSED No CoA moiety '
               'found (missing pantetheine)\n'
               ' * SMILES: '
               '[H]C(CCC(O)=O)=C([H])C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 2,3-didehydroadipoyl-CoA REASON: MISSED No CoA moiety '
               'found (missing pantetheine)\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-2-undecenoyl-CoA REASON: MISSED No CoA moiety '
               'found (missing pantetheine)\n'
               ' * SMILES: '
               'CC=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: but-2-enoyl-CoA REASON: MISSED No CoA moiety found '
               '(missing pantetheine)\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z)-3-oxoicosa-11-enoyl-CoA REASON: MISSED No CoA '
               'moiety found (missing pantetheine)\n'
               ' * SMILES: '
               'CC(C)CCCC(C)CCCC(C)CC\\C=C(/C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (E)-2,3-didehydropristanoyl-CoA REASON: MISSED No CoA '
               'moiety found (missing pantetheine)\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(CCCCC/C=C\\CCCCCCCC)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (7Z)-hexadecenoyl-CoA REASON: MISSED No CoA moiety found '
               '(missing pantetheine)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze what went wrong and suggest improvements:\n'
                 '\n'
                 '1. The main issue seems to be with the CoA moiety detection. '
                 'The error messages show that the code is failing to detect '
                 'the pantetheine part of CoA in valid molecules.\n'
                 '\n'
                 '2. Looking at the example SMILES, I notice that the '
                 'pantetheine pattern we used ("NC(=O)CCNC(=O)CCSC(=O)") might '
                 'be too rigid. The CoA structure can have variations in how '
                 "it's represented in SMILES.\n"
                 '\n'
                 '3. We should also consider that some SMILES strings might '
                 'represent the same structure but with different bond '
                 'orientations (E/Z or cis/trans).\n'
                 '\n'
                 "Here's my improved approach:\n"
                 '\n'
                 '1. Make the CoA pattern detection more flexible by:\n'
                 '   - Breaking it into smaller essential parts\n'
                 '   - Using more general SMARTS patterns\n'
                 '   - Adding phosphate group detection\n'
                 '\n'
                 '2. Improve double bond detection by:\n'
                 '   - Better handling of E/Z notation\n'
                 '   - More precise aromatic vs non-aromatic distinction\n'
                 '\n'
                 "3. Add molecular weight and size checks to ensure we're "
                 'dealing with fatty acid-sized molecules\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 83,
    'num_false_positives': 100,
    'num_true_negatives': 52147,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.453551912568306,
    'recall': 1.0,
    'f1': 0.6240601503759399,
    'accuracy': 0.9980890502579782,
    'negative_predictive_value': 1.0}