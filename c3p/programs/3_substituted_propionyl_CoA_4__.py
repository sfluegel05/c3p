"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 4 negative charges
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = mol.GetSubstructMatches(charge_pattern)
    if len(charge_matches) != 4:
        return False, f"Found {len(charge_matches)} negative charges, need exactly 4"

    # Check for adenine base (more flexible pattern)
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2[nH0]cnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"

    # Check for ribose (simplified pattern)
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(CO)OC1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose moiety found"

    # Check for phosphate groups (more flexible pattern)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups found: {phosphate_matches}"

    # Check for pantetheine arm (more flexible pattern)
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine arm found"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Verify presence of acyl chain with at least 3 carbons
    acyl_chain = Chem.MolFromSmarts("SC(=O)CC")
    if not mol.HasSubstructMatch(acyl_chain):
        return False, "Acyl chain too short or missing"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # Basic sanity checks on atom counts
    if p_count != 3:
        return False, f"Must have exactly 3 phosphorus atoms, found {p_count}"
    if s_count != 1:
        return False, f"Must have exactly 1 sulfur atom, found {s_count}"
    if n_count < 5:
        return False, f"Insufficient nitrogen atoms for CoA structure, found {n_count}"
    if c_count < 20:
        return False, f"Insufficient carbon atoms for CoA structure, found {c_count}"

    # All criteria met
    return True, "Contains CoA(4-) moiety with appropriate thioester linkage and substituted acyl chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65111',
                          'name': '3-substituted propionyl-CoA(4-)',
                          'definition': 'An acyl-CoA(4-) oxoanion arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        '3-substituted propionyl-CoA; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:58342'],
                          'xrefs': ['PMID:17719544'],
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
               'CC(C)CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: isopentadecanoyl-CoA(4-) REASON: MISSED No '
               'ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z,14Z,17Z,20Z,23Z)-hexacosapentaenoyl-CoA(4-) REASON: '
               'MISSED No ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z,14Z,17Z)-icosatrienoyl-CoA(4-) REASON: MISSED No '
               'ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (13Z,16Z,19Z,22Z)-octacosatetraenoyl-CoA(4-) REASON: '
               'MISSED No ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (19Z,22Z,25Z,28Z)-tetratriacontatetraenoyl-CoA(4-) '
               'REASON: MISSED No ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: hexacosanoyl-CoA(4-) REASON: MISSED No ribose-phosphate '
               'moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z,14Z)-icosadienoyl-CoA(4-) REASON: MISSED No '
               'ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: (11Z,14Z,17Z,20Z)-hexacosatetraenoyl-CoA(4-) REASON: '
               'MISSED No ribose-phosphate moiety found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: montanoyl-CoA(4-) REASON: MISSED No ribose-phosphate '
               'moiety found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12 '
               'NAME: icosanoyl-CoA(4-) REASON: MISSED No ribose-phosphate '
               'moiety found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Found 0 negative charges, need '
                                               'exactly 4'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Let me analyze the issues and propose improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the SMARTS patterns, '
                 'especially for the ribose-phosphate moiety, are too specific '
                 'and rigid. The previous code failed to match any of the '
                 'positive examples, indicating the patterns need to be more '
                 'flexible.\n'
                 '\n'
                 '2. The error messages show that all molecules failed at the '
                 '"No ribose-phosphate moiety found" check, suggesting this '
                 'pattern needs to be revised.\n'
                 '\n'
                 '3. Looking at the example structures, they all share common '
                 'features:\n'
                 '   - A CoA backbone with 4 negative charges\n'
                 '   - An adenine nucleobase\n'
                 '   - A pantetheine arm\n'
                 '   - A thioester linkage\n'
                 '   - A variable acyl chain\n'
                 '\n'
                 '4. Improvements needed:\n'
                 '   - Simplify the ribose-phosphate pattern to be less '
                 'restrictive\n'
                 '   - Break down complex patterns into smaller, more flexible '
                 'components\n'
                 '   - Focus on the essential structural elements rather than '
                 'exact matches\n'
                 '   - Add checks for the 3-substituted propionyl group '
                 'specifically\n'
                 "   - Consider stereochemistry but don't make it mandatory in "
                 'patterns\n'
                 '\n'
                 "Here's the improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 30,
    'num_false_positives': 100,
    'num_true_negatives': 35263,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.23076923076923078,
    'recall': 1.0,
    'f1': 0.375,
    'accuracy': 0.9971745825445709,
    'negative_predictive_value': 1.0}