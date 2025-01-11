"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:2468 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone with a hydroxyl group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alpha-hydroxy ketone
    # The pattern looks for a carbonyl group (C=O) with a hydroxyl group (OH) on the alpha-carbon
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][OX2H]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains a ketone with a hydroxyl group on the alpha-carbon"
    else:
        return False, "No alpha-hydroxy ketone pattern found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139588',
                          'name': 'alpha-hydroxy ketone',
                          'definition': 'A ketone containing a hydroxy group '
                                        'on the alpha-carbon relative to the '
                                        'C=O group.',
                          'parents': ['CHEBI:17087', 'CHEBI:33822'],
                          'xrefs': [   'PMID:15326516',
                                       'PMID:19908854',
                                       'PMID:20382022',
                                       'PMID:23295224'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC',
                                     'name': 'TG(22:1(13Z)/16:1(9Z)/18:4(6Z,9Z,12Z,15Z))',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'O1C2C(C3(C(C4C(C5(C(CC4)CC(OC6OC(C(OC7OC(C(O)C(O)C7O)C)C(O)C6O)COC8OC(C(O)C(O)C8O)CO)CC5)C)CC3)C2)C)C(C19OCC(CC9)C)C',
                                     'name': 'Desglucoparillin',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': '[C@H]1(O)C(O)O[C@H](C=O)[C@H]([C@@H]1O)O',
                                     'name': '6-dehydro-D-glucose',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'C1=CC(=CC=C1NC2=C(C=NC=C2)S(=O)(=O)N)Cl',
                                     'name': '4-(4-chloroanilino)-3-pyridinesulfonamide',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'C=1([NH+]=C(C=C(N1)NCC=2ON=C(N2)C(N)=O)C(C)C)N',
                                     'name': '2-amino-4-{[(3-carbamoyl-1,2,4-oxadiazol-5-yl)methyl]amino}-6-isopropylpyrimidin-1-ium',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'O([C@@H]1[C@@H](O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@H](O[C@@H]([C@@H]1O)CO)O)[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O)CO',
                                     'name': '(2R,3S,4S,5S,6R)-2-[(2S,3R,4S,5S,6R)-2,5-Dihydroxy-6-(hydroxymethyl)-3-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-4-yl]oxy-6-(hydroxymethyl)oxane-3,4,5-triol',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)S(=O)(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@@H](O[C@@H]4CO)CC(=O)O',
                                     'name': '2-[(1R,3R,4aS,9aR)-1-(hydroxymethyl)-6-[(4-methoxyphenyl)sulfonylamino]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]acetic '
                                             'acid',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'O=C(C=1NC=CC1)[C@@H]2[C@@H](C=C[C@H]3[C@H]2CCC3)C=CC=C([C@@H]4O[C@@H]([C@H](C(=O)O)C)CC[C@@H]4C)CC',
                                     'name': 'A83094A',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'COc1ccc(\\C=C\\C(=O)N[C@H](CCCNC(N)=N)C(O)=O)cc1Br',
                                     'name': 'subereamine A',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'},
                                 {   'smiles': 'O=C1NCCCCC1NC(=O)C(=O)NCC=2C=CC=NC2',
                                     'name': 'N1-(2-oxoazepan-3-yl)-N2-(3-pyridylmethyl)ethanediamide',
                                     'reason': 'No alpha-hydroxy ketone '
                                               'pattern found'}],
    'sample_false_negatives': [   {   'smiles': 'OC1=CC=CC(O)=C(O)C1=O',
                                      'name': '3,7-dihydroxytropolone',
                                      'reason': 'No alpha-hydroxy ketone '
                                                'pattern found'},
                                  {   'smiles': 'OC1=CC=CC=CC1=O',
                                      'name': 'tropolone',
                                      'reason': 'No alpha-hydroxy ketone '
                                                'pattern found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as an alpha-hydroxy ketone, we need '
                 'to identify the presence of a ketone group (C=O) and a '
                 'hydroxyl group (OH) on the alpha-carbon relative to the '
                 'carbonyl group. The alpha-carbon is the carbon atom directly '
                 'adjacent to the carbonyl carbon.\n'
                 '\n'
                 'Here’s the step-by-step reasoning:\n'
                 '1. **Identify the ketone group**: We look for a carbonyl '
                 'group (C=O) in the molecule.\n'
                 '2. **Identify the alpha-carbon**: The alpha-carbon is the '
                 'carbon directly adjacent to the carbonyl carbon.\n'
                 '3. **Check for a hydroxyl group on the alpha-carbon**: We '
                 'verify if the alpha-carbon has a hydroxyl group (OH) '
                 'attached to it.\n'
                 '\n'
                 'The function `is_alpha_hydroxy_ketone` will implement this '
                 'logic using RDKit.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 405,
    'num_false_positives': 100,
    'num_true_negatives': 1783,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.801980198019802,
    'recall': 0.995085995085995,
    'f1': 0.888157894736842,
    'accuracy': 0.9554585152838428,
    'negative_predictive_value': 0.9988795518207283}