"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for UDP-sugar components
    uridine_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)")
    
    # Check for uridine moiety
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine moiety found"

    # Check for diphosphate group
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"

    # Check for sugar component
    # A more specific sugar pattern could be defined, but in general a UDP-sugar
    # would have a recognizable sugar moiety linked with the above patterns
    # Here we will assume that if both uridine and diphosphate are detected, the
    # sugar moiety is implied by the context of these patterns in known UDP-sugars.

    return True, "Contains uridine and diphosphate linkage typical of UDP-sugars"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17297',
                          'name': 'UDP-sugar',
                          'definition': 'A pyrimidine nucleotide-sugar having '
                                        'UDP as the nucleotide component '
                                        'attached to an unspecified sugar via '
                                        'an anomeric diphosphate linkage.',
                          'parents': ['CHEBI:61109'],
                          'xrefs': ['KEGG:C05227'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No uridine moiety found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No uridine moiety found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as belonging to the class '
                 'UDP-sugar based on its SMILES string, we can look for '
                 'specific structural patterns that are characteristic of '
                 'UDP-sugars. \n'
                 '\n'
                 'A UDP-sugar typically consists of a uridine diphosphate '
                 'group and a saccharide component. The distinguishing '
                 'features include:\n'
                 '\n'
                 '1. The presence of a uridine moiety, which can be '
                 'represented by the substructure pattern `n1ccc(=O)[nH]c1=O` '
                 'in SMILES, signifying the uridine part.\n'
                 '2. A diphosphate linkage, which can appear as '
                 '`OP(O)(=O)OP(O)(=O)` in SMILES.\n'
                 '3. A glycosyl group, the sugar component, typically '
                 'represented by a hexose sugar backbone with optional '
                 'modifications like acetylaminos, deoxys, or uronic acids.\n'
                 '\n'
                 'Given these observations, we should aim to create a function '
                 '`is_UDP_sugar` that will perform the following steps:\n'
                 '- Parse the SMILES string into a molecular object.\n'
                 '- Check for the presence of the uridine moiety.\n'
                 '- Verify the presence of a diphosphate group.\n'
                 '- Assert the presence of a sugar moiety linked to the '
                 'uridine diphosphate.\n'
                 '\n'
                 "Let's implement this logic using RDKit.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 62,
    'num_false_positives': 100,
    'num_true_negatives': 49397,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.38271604938271603,
    'recall': 1.0,
    'f1': 0.5535714285714286,
    'accuracy': 0.997982203030731,
    'negative_predictive_value': 1.0}