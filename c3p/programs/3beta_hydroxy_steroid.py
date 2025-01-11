"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid
A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for 3-OH group with specific stereochemistry
    # [OH] connected to carbon at position 3 with beta stereochemistry
    # The [H] indicates explicit hydrogen, helps define stereochemistry
    # The @@ indicates beta configuration (below the plane)
    beta_3_oh_pattern = Chem.MolFromSmarts('[H][C@@]1[C@@H](O)CC[C@]2')
    
    # Alternative pattern that might match other valid 3beta-OH steroids
    alt_beta_3_oh_pattern = Chem.MolFromSmarts('[C@@H](O)CC[C@@]1')
    
    if not (mol.HasSubstructMatch(beta_3_oh_pattern) or mol.HasSubstructMatch(alt_beta_3_oh_pattern)):
        return False, "No 3beta-hydroxy group found"

    # Count rings to ensure we have a steroid-like structure
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check carbon count (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    return True, "Contains steroid core with 3beta-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36836',
                          'name': '3beta-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the '
                                        '3-hydroxy substituent is in the  '
                                        'beta-position.',
                          'parents': ['CHEBI:35681', 'CHEBI:36834'],
                          'xrefs': [   'KEGG:C02945',
                                       'MetaCyc:3-Beta-Hydroxysterols',
                                       'PMID:10535978',
                                       'PMID:12829805'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCC(=C)C(C)C',
                                      'name': 'ergosta-5,7,24(28)-trien-3beta-ol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])C(C)(C)[C@@H](O)CC[C@@]44C[C@@]34CC[C@]12C)[C@H](C)CC[C@@H](C)C(C)=C',
                                      'name': '(24R)-24-methylcycloart-25-en-3beta-ol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)[C@@H](C=O)[C@@H]1CC3',
                                      'name': '4alpha-formyl-5alpha-cholest-8-en-3beta-ol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@H](O)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CC[C@H](O)C2',
                                      'name': '5beta-androstane-3beta,17beta-diol',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC(C)C(=C)CC[C@@H](C)[C@H]1CC[C@H]2C3=C([C@H](O)C[C@]12C)[C@@]1(C)CC[C@H](O)C[C@@H]1CC3=O',
                                      'name': '3beta,11alpha-dihydroxyergosta-8,24(28)-dien-7-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12C[C@@]3([H])[C@]4([H])CC=C5C[C@@H](O)CC[C@]5(C)[C@@]4([H])CC[C@]3(C)[C@@]1([H])[C@H](C)[C@@]1(CC[C@@](C)(CO)O1)O2',
                                      'name': 'nuatigenin',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@]3([H])[C@]([H])(C[C@@H](O)[C@]4(C)[C@]([H])(CC[C@]34O)C3=CC(=O)OC3)[C@@]1(C)CC[C@H](O)C2',
                                      'name': 'digoxigenin',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)[C@H](O)C[C@H]1C[C@@H](O)CC[C@]31C',
                                      'name': '(3beta,5alpha,7alpha)-3,7-dihydroxyandrostan-17-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[C@]12([C@]3([C@]([H])([C@@H](C[C@]2([H])[C@@]4([C@@]([H])([C@H](C1)O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)CO)C([C@H](CC4)O)(C)C)C)O)[C@]([H])(CC3)[C@](CCC=C(C)C)(O)C)C)C',
                                      'name': '(20S)-ginsenoside Rh1',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'CC(C)[C@@H](C)\\C=C\\[C@@H](C)[C@H]1CC[C@H]2C3=CC(=O)[C@@]4(O)C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': '(22E,24R)-3beta,5alpha-dihydroxyergosta-7,22-dien-6-one',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify 3beta-hydroxy steroids. "
                 'Looking at the examples, these compounds share several key '
                 'characteristics:\n'
                 '\n'
                 '1. They have a steroid core structure\n'
                 '2. They have a hydroxyl group (-OH) at the 3-position\n'
                 '3. The 3-OH is in the beta configuration (below the plane of '
                 'the steroid)\n'
                 '\n'
                 "Here's the program:",
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