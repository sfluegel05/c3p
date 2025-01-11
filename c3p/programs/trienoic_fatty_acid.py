"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for three or more double bonds; consider cis/trans isomerism
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found {len(double_bond_matches)} double bonds, need at least 3 for trienoic structure"
    
    # Count number of carbon atoms to check for long chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Typically, a trienoic fatty acid has a moderate to long chain
    if c_count < 12:
        return False, "Too short hydrocarbon chain for a typical trienoic fatty acid"
    
    return True, "Contains carboxylic acid group and at least three double bonds, typical of trienoic fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73155',
                          'name': 'trienoic fatty acid',
                          'definition': 'Any polyunsaturated fatty acid that '
                                        'contains three double bonds.',
                          'parents': ['CHEBI:26208', 'CHEBI:53339'],
                          'xrefs': ['PMID:11115886'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found 1 double bonds, need at '
                                               'least 3 for trienoic '
                                               'structure'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Found 0 double bonds, need at '
                                               'least 3 for trienoic '
                                               'structure'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found 0 double bonds, need at '
                                               'least 3 for trienoic '
                                               'structure'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Found 0 double bonds, need at '
                                               'least 3 for trienoic '
                                               'structure'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'Found 0 double bonds, need at '
                                               'least 3 for trienoic '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1N(C(=O)C=2N=CN=NC2N1C)C',
                                      'name': 'Fervenulin',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(CCC1(CC)C(=O)NC(=O)NC1=O)(C)C',
                                      'name': "3'-Hydroxyamobarbital",
                                      'reason': 'No carboxylic acid group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as trienoic fatty acids based '
                 'on their SMILES strings, we need to identify key structural '
                 'features characteristic of this class. A trienoic fatty acid '
                 'typically has:\n'
                 '\n'
                 '1. A carboxylic acid group (`-COOH`) at one end of the '
                 'molecule.\n'
                 '2. At least three carbon-carbon double bonds, which are '
                 'usually conjugated and can have various stereochemical '
                 'configurations (cis/trans).\n'
                 '3. Typically a long hydrocarbon chain resembling that of '
                 'fatty acids, generally linear but may have some functional '
                 'groups or rings.\n'
                 '\n'
                 'Given these requirements, we will implement a function '
                 '`is_trienoic_fatty_acid` that will check for these features '
                 'in a SMILES string using RDKit, a cheminformatics library.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 135,
    'num_false_positives': 100,
    'num_true_negatives': 1469,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.574468085106383,
    'recall': 0.9854014598540146,
    'f1': 0.7258064516129032,
    'accuracy': 0.94021101992966,
    'negative_predictive_value': 0.9986403806934059}