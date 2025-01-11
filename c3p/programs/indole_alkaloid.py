"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton and typically has nitrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded indole patterns to match more variations
    indole_patterns = [
        Chem.MolFromSmarts("[nX3]1ccc2ccccc12"),  # Basic indole with any nitrogen type
        Chem.MolFromSmarts("[nX2]1ccc2ccccc12"),  # Indole with double-bonded nitrogen
        Chem.MolFromSmarts("[n+]1ccc2ccccc12"),   # Indole with charged nitrogen
        Chem.MolFromSmarts("[nH]1ccc2ccccc12"),   # Indole with hydrogen on nitrogen
        Chem.MolFromSmarts("[nX3]1ccc2c1cccc2"),  # Indole with fused rings
        Chem.MolFromSmarts("[nX2]1ccc2c1cccc2"),  # Indole with fused rings and double-bonded nitrogen
        Chem.MolFromSmarts("[n+]1ccc2c1cccc2"),   # Indole with fused rings and charged nitrogen
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2")    # Indole with fused rings and hydrogen on nitrogen
    ]

    # Check for any indole pattern match
    indole_found = any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns)
    if not indole_found:
        return False, "No indole skeleton found"

    # Check for nitrogen atoms in the context of the indole skeleton
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count == 0:
        return False, "No nitrogen atoms found"

    # Check if at least one nitrogen is part of the indole skeleton
    indole_nitrogen_found = any(
        any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7 for atom_idx in match)
        for pattern in indole_patterns
        for match in mol.GetSubstructMatches(pattern)
    )
    if not indole_nitrogen_found:
        return False, "No nitrogen atoms in indole skeleton"

    # Check molecular weight - indole alkaloids typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for indole alkaloid"

    return True, "Contains indole skeleton and nitrogen atoms"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38958',
                          'name': 'indole alkaloid',
                          'definition': 'An alkaloid containing an indole skeleton.',
                          'parents': ['CHEBI:38958', 'CHEBI:22315']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38958',
                          'name': 'indole alkaloid',
                          'definition': 'An alkaloid containing an indole '
                                        'skeleton.',
                          'parents': ['CHEBI:22315'],
                          'xrefs': ['KEGG:C06073', 'Wikipedia:Indole_alkaloid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: '<' not supported between instances of 'slice' and "
               "'int'\n"
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12C[C@]3([H])N([C@H]4C[C@@]5([C@H](OC(C)=O)C14)c1ccccc1N[C@@]35[H])[C@H](O)\\C2=C\\C',
                                      'name': '1,2-dihydrovomilenine',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'C/C=C\\1/CN2CC[C@]34C5=CC=CC=C5N=C4[C@@]2(C[C@@]1(C3(COC(C)=O)C(=O)OC)[H])[H]',
                                      'name': 'Akuammiline',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': '[C@@]123[C@@](N(C4=C1C=CC(=C4)OC)C)([C@]([C@@H]([C@]5([C@@]2(N(CC=C5)CC3)[H])CC)OC(=O)C)(C(=O)OC)O)[H]',
                                      'name': 'vindoline',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': '[H][C@@]12Nc3ccccc3[C@@]1(CCN2C)[C@]12CCN(C)[C@@]1([H])Nc1ccccc21',
                                      'name': 'meso-chimonanthine',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'O[C@H]1[C@]23[C@]([C@]4([N@@+]5([C@@H]([C@H]([C@@H]([C@H]1[C@]5([H])C2)C4)CC)O)CC(CN(CC)CC)O)[H])(N(C=6C=CC=CC36)C)[H]',
                                      'name': 'detajmium',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'C=12[C@@]3([C@H]4[NH+](CC3)C\\C(\\[C@@](C1C=O)(C4)[H])=C\\C)C5=CC=CC=C5N2',
                                      'name': 'norfluorocurarine',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'COC([C@@H]1C[C@@]23CCCN4CCC5(C6=CC=CC=C6N(C)[C@]15CC2)[C@]34[H])=O',
                                      'name': 'Pleiocarpinine',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'N12C([C@@H](NC([C@@]1(CCC2)[H])=O)C[C@]3(C(NC4=C3C=CC5=C4C=CC(O5)(C)C)=O)C(C=C)(C)C)=O',
                                      'name': 'notoamide C',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'C1=CC=C2C(=C1)[C@]34CCN5CC6=CCO[C@H]([C@]([C@]6(C[C@@]35[H])[H])([C@@]4(N2)[H])[H])O',
                                      'name': 'caracurine VII',
                                      'reason': 'No indole skeleton found'},
                                  {   'smiles': 'C/C=C\\1/CN2[C@](C[C@@]34C5=CC=CC=C5N6C(C=C[C@@]([C@]1(C[C@@]32[H])[H])([C@@]46[H])[H])=O)([C@@]7(C=C8[C@]/9(C[C@]%10([C@@]%11(CCN%10C\\C9=C\\CO)C%12=CC=CC=C%12N([C@@]8%11[H])C7=O)[H])[H])[H])[H]',
                                      'name': '18-hydroxyisosungucine',
                                      'reason': 'No indole skeleton found'}],
    'attempt': 3,
    'reasoning': 'The previous code failed because of an error in the way the '
                 'nitrogen atoms in the indole skeleton were being checked. '
                 'Specifically, the line `for atom in '
                 'mol.GetAtoms()[match[0]:match[-1]+1]` caused a TypeError '
                 'because `match` is a tuple of atom indices, and slicing a '
                 'tuple with a slice object is not supported.\n'
                 '\n'
                 'To fix this, we need to iterate over the atom indices in the '
                 'match and check if any of them correspond to a nitrogen '
                 'atom. Additionally, we can simplify the logic by directly '
                 'checking if the matched substructure contains a nitrogen '
                 'atom.\n'
                 '\n'
                 'Here is the corrected code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 180,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.6428571428571429,
    'f1': 0.782608695652174,
    'accuracy': 0.6428571428571429,
    'negative_predictive_value': 0.0}