"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety
    # Look for adenine base
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety of CoA"
    
    # Look for ribose phosphate
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c2ncnc3)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose phosphate portion of CoA"
    
    # Look for pantetheine portion with thioester
    # [CX4]C(=O)NCCC(=O)NCCS represents the pantetheine chain
    # The C(=O)S represents the thioester linkage
    pantetheine_thioester = Chem.MolFromSmarts("[CX4]C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_thioester):
        return False, "Missing pantetheine-thioester portion"
    
    # Check for 3-oxo group pattern
    # Pattern: carbon chain - C(=O)CC(=O)S- 
    # This represents the 3-oxo group followed by the thioester
    oxo_pattern = Chem.MolFromSmarts("C-C(=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group pattern"
    
    # Additional checks for fatty acid portion
    # Count carbons in the main chain
    # We'll be lenient here as fatty acids can vary in length
    carbon_chain = Chem.MolFromSmarts("CCCCC")  # At least 5 carbons
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Carbon chain too short for fatty acid portion"

    # Check for reasonable molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 800:  # CoA itself is quite large
        return False, "Molecular weight too low for 3-oxo-fatty acyl-CoA"
    
    return True, "Contains CoA moiety, 3-oxo group, and appropriate fatty acid chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15489',
                          'name': '3-oxo-fatty acyl-CoA',
                          'definition': 'An oxo fatty acyl-CoA that results '
                                        'from the formal condensation of the '
                                        'thiol group of coenzyme A with the '
                                        'carboxy group of any 3-oxo-fatty '
                                        'acid..',
                          'parents': ['CHEBI:61903'],
                          'xrefs': [   'KEGG:C00264',
                                       'PMID:11315193',
                                       'PMID:11418601',
                                       'PMID:11879205',
                                       'PMID:7957058',
                                       'PMID:8541311'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'Missing adenine moiety of CoA'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'Missing adenine moiety of '
                                               'CoA'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help you create a program to classify 3-oxo-fatty '
                 'acyl-CoA compounds. Looking at the examples, these compounds '
                 'have several key characteristics:\n'
                 '\n'
                 '1. They contain a CoA (Coenzyme A) moiety\n'
                 '2. They have a 3-oxo group (ketone at position 3)\n'
                 '3. They have a thioester linkage connecting the fatty acid '
                 'to CoA\n'
                 '4. The fatty acid portion can be saturated or unsaturated, '
                 'and varies in length\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 59,
    'num_false_positives': 100,
    'num_true_negatives': 95099,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3710691823899371,
    'recall': 1.0,
    'f1': 0.5412844036697247,
    'accuracy': 0.9989502194041445,
    'negative_predictive_value': 1.0}