"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing a carbonyl group bonded to one hydrogen atom and one R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the aldehyde group: R-C(=O)H
    aldehyde_pattern = Chem.MolFromSmarts("[#6]=O")

    # Check for the aldehyde group in the molecule
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde group (R-C(=O)H)"
    else:
        return False, "No aldehyde group found"

# Example usage
example_smiles = "C=O"  # Formyl group in SMILES
result, reason = is_aldehyde(example_smiles)
print(f"Is aldehyde: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17478',
                          'name': 'aldehyde',
                          'definition': 'A compound RC(=O)H, in which a '
                                        'carbonyl group is bonded to one '
                                        'hydrogen atom and to one R group.',
                          'parents': ['CHEBI:36586'],
                          'xrefs': ['KEGG:C00071'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'O1C([C@@H](O)C[C@@H]1[C@H]([C@@H]2[C@@]3([C@@](C=4C([C@@]5(C(C([C@@H](O)CC5)(C)C)CC4)C)=CC3)(C)CC2)C)C)(C)C',
                                     'name': 'Ganoderchochlearin A',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'N=1C2=C(N(C=3C=CC=CC3C2=CC1)C)C',
                                     'name': 'Marinoquinoline G',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': 'Clc1ccc(CN2CCS\\C2=N/C#N)cn1',
                                     'name': '(Z)-thiacloprid',
                                     'reason': 'No aldehyde group found'},
                                 {   'smiles': '[Na+].[Na+].[O-]P(=O)(O)[O-].O.O.O.O.O.O.O.O.O.O.O.O',
                                     'name': 'disodium hydrogenphosphate '
                                             'dodecahydrate',
                                     'reason': 'No aldehyde group found'}],
    'sample_false_negatives': [   {   'smiles': 'C1C[C@H](N(C1)C(=[NH2+])N)O',
                                      'name': '(2R)-2-hydroxypyrrolidine-1-carboximidamide(1+)',
                                      'reason': 'No aldehyde group found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as an aldehyde, we need to check for '
                 'the presence of a carbonyl group (C=O) bonded to a hydrogen '
                 'atom and to one R group. We can use SMARTS patterns in RDKit '
                 'to identify this substructure.\n'
                 '\n'
                 'The general strategy will involve:\n'
                 '\n'
                 '1. Parsing the SMILES string into a molecule object using '
                 'RDKit.\n'
                 '2. Defining a SMARTS pattern that matches the aldehyde '
                 'functional group, RC(=O)H.\n'
                 '3. Using the `HasSubstructMatch` method in RDKit to check if '
                 'the aldehyde pattern exists within the molecule.\n'
                 '4. Returning `True` with a reason if the pattern matches, '
                 'otherwise returning `False` with a reason.\n'
                 '\n'
                 "Here's the Python program that accomplishes the task:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 787,
    'num_false_positives': 100,
    'num_true_negatives': 19,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.887260428410372,
    'recall': 0.998730964467005,
    'f1': 0.9397014925373134,
    'accuracy': 0.8886438809261301,
    'negative_predictive_value': 0.95}