"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ketone group pattern: a carbonyl group (C=O) with carbon atoms on either side
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Find substructure matches for ketone groups
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check for at least two ketone groups
    if len(ketone_matches) >= 2:
        return True, f"Contains {len(ketone_matches)} ketone groups"

    return False, f"Found {len(ketone_matches)} ketone groups, need at least 2 for diketone"

# Example usage:
# print(is_diketone("O=C(CCCCCCCCCCCCCCCCC)CC(=O)CCCCCC"))  # Hexacosane-7,9-dione


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46640',
                          'name': 'diketone',
                          'definition': 'A compound that contains two ketone '
                                        'functionalities.',
                          'parents': ['CHEBI:17087'],
                          'xrefs': ['Wikipedia:Diketone'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'Found 1 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'Found 1 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine hydrochloride',
                                     'reason': 'Found 0 ketone groups, need at '
                                               'least 2 for diketone'},
                                 {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=CC(=C3)C#N)O[C@H]1CN(C)C(=O)C4=CC=NC=C4)[C@H](C)CO',
                                     'name': 'N-[[(4R,5R)-8-(3-cyanophenyl)-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-N-methyl-4-pyridinecarboxamide',
                                     'reason': 'Found 1 ketone groups, need at '
                                               'least 2 for diketone'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as a diketone, we need to look '
                 'for a specific structural motif: two ketone groups within '
                 'the molecule. A ketone group is characterized by a carbonyl '
                 'group (C=O) bonded to two alkyl groups. To simplify the '
                 'task, we can use the RDKit library, which allows us to work '
                 'with SMILES strings and search for substructures using '
                 'SMARTS patterns. The following steps outline the approach:\n'
                 '\n'
                 '1. Parse the SMILES string to create a molecular object '
                 "using RDKit's `MolFromSmiles` function.\n"
                 '2. Define a SMARTS pattern for a ketone group. The pattern '
                 'can be defined as `[CX3](=O)[#6]`, which represents a carbon '
                 'atom double bonded to an oxygen (carbonyl group) and single '
                 'bonded to another carbon atom.\n'
                 "3. Use RDKit's `GetSubstructMatches` method to find matches "
                 'of the defined ketone pattern in the molecule.\n'
                 '4. Verify if the molecule contains at least two distinct '
                 'ketone groups by checking the number of matches.\n'
                 '5. Return `True` and a reason if the molecule contains two '
                 'or more ketone groups, otherwise return `False` with an '
                 'appropriate reason.\n'
                 '\n'
                 'Here is the Python code for the `is_diketone` function:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 155,
    'num_false_positives': 100,
    'num_true_negatives': 145,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6078431372549019,
    'recall': 1.0,
    'f1': 0.7560975609756097,
    'accuracy': 0.75,
    'negative_predictive_value': 1.0}