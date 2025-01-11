"""
Classifies: CHEBI:60834 beta-carbolines
"""
"""
Classifies compounds as beta-carbolines based on their SMILES structure
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines contain a pyrido[3,4-b]indole core structure and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for beta-carboline core structure
    # Pattern 1: Fully aromatic beta-carboline core
    beta_carboline_pattern1 = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2nccc3")
    
    # Pattern 2: Partially hydrogenated beta-carboline core (allowing for variations in saturation)
    beta_carboline_pattern2 = Chem.MolFromSmarts("C1CCc2c(C1)[nH]c3c2nccc3")
    beta_carboline_pattern3 = Chem.MolFromSmarts("[#6]1[#6][#6]c2c([#6]1)[nH]c3c2[#7][#6][#6][#6]3")
    
    # Pattern 4: N-substituted variants (common in the examples)
    beta_carboline_pattern4 = Chem.MolFromSmarts("c1ccc2c(c1)n([#6,#1])c3c2nccc3")
    beta_carboline_pattern5 = Chem.MolFromSmarts("[#6]1[#6][#6]c2c([#6]1)n([#6,#1])c3c2[#7][#6][#6][#6]3")

    # Check for any of the patterns
    if mol.HasSubstructMatch(beta_carboline_pattern1):
        return True, "Contains fully aromatic beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern2):
        return True, "Contains partially hydrogenated beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern3):
        return True, "Contains modified beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern4):
        return True, "Contains N-substituted beta-carboline core structure"
    elif mol.HasSubstructMatch(beta_carboline_pattern5):
        return True, "Contains modified N-substituted beta-carboline core structure"

    # Additional checks for identifying potential false negatives
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient ring systems for beta-carboline structure"

    return False, "Does not contain beta-carboline core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60834',
                          'name': 'beta-carbolines',
                          'definition': 'Any pyridoindole containing a '
                                        'beta-carboline skeleton and their '
                                        'hydrogenated derivatives',
                          'parents': ['CHEBI:48888'],
                          'xrefs': ['PMID:1967844'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](N(CC34CN(C4)C(=O)C5=CC=CC=C5F)CC6=CC=CC=C6)CO',
                                      'name': "(2-fluorophenyl)-[(1S)-1-(hydroxymethyl)-7-methoxy-9-methyl-2-(phenylmethyl)-1'-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanone",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](N(CC34CCN(CC4)CC5CCCC5)CC6=NC=CS6)CO',
                                      'name': "[(1S)-1'-(cyclopentylmethyl)-7-methoxy-9-methyl-2-(2-thiazolylmethyl)-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CC1=C(C(=NO1)C)NC(=O)N2CCC3(CC2)CN([C@H](C4=C3C5=C(N4)C=C(C=C5)OC)CO)CC6=CC7=C(C=C6)OCO7',
                                      'name': "(1R)-2-(1,3-benzodioxol-5-ylmethyl)-N-(3,5-dimethyl-4-isoxazolyl)-1-(hydroxymethyl)-7-methoxy-1'-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]carboxamide",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N2CC3(C2)CN([C@H](C4=C3C5=C(N4C)C=C(C=C5)OC)CO)CC6=CC=CC=C6',
                                      'name': "[(1R)-7-methoxy-9-methyl-1'-(4-methylphenyl)sulfonyl-2-(phenylmethyl)-1-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]methanol",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CCS(=O)(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3)C=C(C=C4)OC)CO)C(=O)NC5=C(ON=C5C)C',
                                      'name': "(1R)-N-(3,5-dimethyl-4-isoxazolyl)-1'-ethylsulfonyl-1-(hydroxymethyl)-7-methoxy-2-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]carboxamide",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'COC1=CC2=C(C=C1)C3=C(N2)[C@H](N(CC34CN(C4)C(=O)C5=CC=CC=C5)C(=O)CC6CC6)CO',
                                      'name': "1-[(1S)-1'-benzoyl-1-(hydroxymethyl)-7-methoxy-2-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,3'-azetidine]yl]-2-cyclopropylethanone",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'COc1c(O)c(=O)n2c3ccccc3c3ccnc1c23',
                                      'name': 'Nigakinone',
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'COC1=CC2=C(C=C1)C3=C(N2)[C@@H](N(CC34CCN(CC4)S(=O)(=O)C5=CC=CS5)C(=O)C6=CC=C(C=C6)F)CO',
                                      'name': "(4-fluorophenyl)-[(1R)-1-(hydroxymethyl)-7-methoxy-1'-thiophen-2-ylsulfonyl-2-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanone",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CN1C2=C(C=CC(=C2)OC)C3=C1[C@@H](NCC34CCN(CC4)CC5=NOC(=C5)C6=CC=CC=C6)CO',
                                      'name': "[(1R)-7-methoxy-9-methyl-1'-[(5-phenyl-3-isoxazolyl)methyl]-1-spiro[2,3-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'},
                                  {   'smiles': 'CC(=O)N1CC2(CN(C2)CC3=CC=C(C=C3)Cl)C4=C([C@H]1CO)N(C5=C4C=CC(=C5)OC)C',
                                      'name': "1-[(1S)-1'-[(4-chlorophenyl)methyl]-1-(hydroxymethyl)-7-methoxy-9-methyl-2-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]yl]ethanone",
                                      'reason': 'Does not contain '
                                                'beta-carboline core '
                                                'structure'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify beta-carbolines. '
                 'Looking at the examples and definition, beta-carbolines are '
                 'characterized by:\n'
                 '\n'
                 '1. A pyrido[3,4-b]indole core structure (also known as '
                 'β-carboline or 9H-pyrido[3,4-b]indole)\n'
                 '2. Can include hydrogenated derivatives\n'
                 '3. The core structure consists of a pyridine ring fused to '
                 'an indole system\n'
                 '\n'
                 "Here's the program to classify beta-carbolines:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0,
    'negative_predictive_value': 0.0}