"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains carbon and at least one chlorine atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of chlorine atoms
    has_chlorine = any(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    
    if has_chlorine:
        return True, "Contains chlorine atoms, classified as an organochlorine compound"
    else:
        return False, "No chlorine atoms found, not an organochlorine compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36683',
                          'name': 'organochlorine compound',
                          'definition': 'An organochlorine compound is a '
                                        'compound containing at least one '
                                        'carbon-chlorine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:23117'],
                          'xrefs': [   'MetaCyc:Chlorides',
                                       'Wikipedia:Organochloride'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'No chlorine atoms found, not '
                                               'an organochlorine compound'}],
    'sample_false_negatives': [   {   'smiles': 'C(C(CO)O)N1C=C(I)C(C(=C1)I)=O',
                                      'name': 'iopydol',
                                      'reason': 'No chlorine atoms found, not '
                                                'an organochlorine compound'},
                                  {   'smiles': 'O=C(N[C@H](C(=O)O)CC1=CC=CC=C1)[C@H]2N(C(=O)[C@@H](NC(=O)[C@@H](N)[C@H](CC)C)CC=3NC=NC3)CCC2',
                                      'name': 'Angiotensin (5-8)',
                                      'reason': 'No chlorine atoms found, not '
                                                'an organochlorine compound'},
                                  {   'smiles': 'O=c1cc(oc2c(cccc12)-c1ccccc1)N1CCOCC1',
                                      'name': 'LY294002',
                                      'reason': 'No chlorine atoms found, not '
                                                'an organochlorine compound'}],
    'attempt': 0,
    'reasoning': 'To create a program that classifies chemical entities as '
                 'organochlorine compounds based on their SMILES strings, we '
                 'need to identify whether a given molecule contains chlorine '
                 'atoms, which is a defining feature of organochlorine '
                 'compounds. Here is the step-by-step process to build the '
                 'program:\n'
                 '\n'
                 '1. **Import necessary libraries**: We will use the RDKit '
                 'library to parse and analyze SMILES strings.\n'
                 '2. **Define the function**: `is_organochlorine_compound` '
                 'will take a SMILES string as input and return a boolean '
                 'value and a reason.\n'
                 '3. **Molecule Parsing**: Use RDKit to convert the SMILES '
                 'string into a molecule object. This will allow us to analyze '
                 'the structure.\n'
                 '4. **Identify Chlorine Atoms**: Check if the molecule '
                 'contains at least one chlorine atom.\n'
                 '5. **Validation and Result**: Return True if the molecule is '
                 'classified as an organochlorine compound with an appropriate '
                 'message; otherwise, return False with a reason.\n'
                 '\n'
                 'Here is the Python code for the function:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1785,
    'num_false_positives': 100,
    'num_true_negatives': 2194,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.946949602122016,
    'recall': 0.9983221476510067,
    'f1': 0.9719575279063435,
    'accuracy': 0.9747672709456149,
    'negative_predictive_value': 0.9986345015930814}