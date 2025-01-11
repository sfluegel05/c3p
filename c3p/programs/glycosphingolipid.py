"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid with a carbohydrate attached via glycosidic linkage 
    to a sphingoid/ceramide base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ceramide/sphingoid base patterns
    # NH-C(=O) amide group
    amide_pattern = Chem.MolFromSmarts("[NH][C;X3](=[O;X1])")
    # Long carbon chain after amide
    chain_pattern = Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found (required for ceramide/sphingoid base)"
    
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found (required for ceramide/sphingoid base)"

    # Look for sugar patterns
    # Pyranose ring with multiple OH groups
    sugar_pattern = Chem.MolFromSmarts("[C;R1][O;R1][C;R1][C;R1][C;R1][C;R1]")
    # Multiple OH groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring found"
    
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, "Not enough hydroxyl groups for sugar moiety"

    # Look for glycosidic linkage (C-O-C)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R1][O;X2][C;X4]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Additional check for sphingoid characteristics
    # Look for OH group near amide and long chain
    sphingoid_oh_pattern = Chem.MolFromSmarts("[NH]C([C;X4])[C;X4][OH]")
    if not mol.HasSubstructMatch(sphingoid_oh_pattern):
        return False, "Missing characteristic sphingoid OH group"

    # Count carbons and oxygens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glycosphingolipid"
    if o_count < 6:
        return False, "Too few oxygens for glycosphingolipid"

    return True, "Contains carbohydrate residue attached via glycosidic linkage to sphingoid/ceramide base"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24402',
                          'name': 'glycosphingolipid',
                          'definition': 'A glycosphingolipid is a glycolipid '
                                        'that is a carbohydrate-containing '
                                        'derivative of a sphingoid or '
                                        'ceramide. It is understood that the '
                                        'carbohydrate residue is attached by a '
                                        'glycosidic linkage to O-1 of the '
                                        'sphingoid.',
                          'parents': ['CHEBI:26739', 'CHEBI:33563'],
                          'xrefs': [   'PMID:18845618',
                                       'PMID:19303901',
                                       'PMID:26954668',
                                       'PMID:26976737',
                                       'PMID:27253729',
                                       'PMID:27325407',
                                       'Wikipedia:Glycosphingolipid'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No long carbon chain found '
                                               '(required for '
                                               'ceramide/sphingoid base)'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No long carbon chain found '
                                               '(required for '
                                               'ceramide/sphingoid base)'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No long carbon chain found '
                                               '(required for '
                                               'ceramide/sphingoid base)'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'No amide group found (required '
                                               'for ceramide/sphingoid base)'}],
    'sample_false_negatives': [   {   'smiles': '[C@H]1([C@H](O[C@@H](OC[C@@H]([C@@H](CCCCCCCCCCCCCCC)O)OC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)[C@@H]([C@H]1O)O)CO)O',
                                      'name': "beta-D-glucosyl-(1<->1')-N-[(15Z)-tetracosenoyl]sphinganine",
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCCC(=S)N[C@@H](CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)[C@H](O)CCCCCCCCCCCCCC',
                                      'name': 'N-[(2S,3S,4R)-1-(alpha-D-galactosyloxy)-3,4-dihydroxyoctadecan-2-yl]hexacosanethioamide',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'N([C@@H](CO[C@H]1[C@@H]([C@@H](OC(=O)C)[C@H]([C@@H](COC(=O)C)O1)OC(=O)C)OC(=O)C)[C@H](OC(=O)C)/C=C/CCCCCCCCCCCCC)C(CCCCCCCCCCCCCCCCCCCCC)=O',
                                      'name': "N-(docosanoyl)-1-(2',3',4',6'-tetraacetylbeta-D-galactosyl)-3-acetylsphingosine",
                                      'reason': 'Not enough hydroxyl groups '
                                                'for sugar moiety'},
                                  {   'smiles': '[H][C@@]1(O[C@H](CCCCCCCCCCCC\\C=C\\C(CCCCCCCC2NCCc3cc(O)c(O)cc23)OC)[C@@H](C)N)O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'oceanalin A',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'psychosine',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CCCCCCCCCC)CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@H]1NC(C)=O',
                                      'name': 'Le(x)-C12L',
                                      'reason': 'Missing characteristic '
                                                'sphingoid OH group'},
                                  {   'smiles': '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O[C@@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O)O)O)O)O)O)OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)[NH3+]',
                                      'name': 'alpha-D-galactosyl-D-galactosyl-D-glucosylsphingosine',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](N)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O',
                                      'name': '1-(3-O-sulfo-beta-D-galactosyl)sphingosine',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)C1OCC(N)C(O)/C=C/CCCCCCCCCCCCC)CO',
                                      'name': '(3R,4S,5S,6R)-2-{[(4E)-2-amino-3-hydroxyoctadec-4-en-1-yl]oxy}-6-(hydroxymethyl)oxane-3,4,5-triol',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'},
                                  {   'smiles': 'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@@H](N)CO[C@@H]1O[C@H](COS(O)(=O)=O)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': 'psychosine sulfate',
                                      'reason': 'No amide group found '
                                                '(required for '
                                                'ceramide/sphingoid base)'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify glycosphingolipids. The "
                 'key characteristics to look for are:\n'
                 '1. A sphingoid/ceramide base\n'
                 '2. A carbohydrate residue\n'
                 '3. A glycosidic linkage to O-1 of the sphingoid\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 335,
    'num_false_positives': 100,
    'num_true_negatives': 5201,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.7701149425287356,
    'recall': 0.9710144927536232,
    'f1': 0.8589743589743589,
    'accuracy': 0.9805171803046404,
    'negative_predictive_value': 0.9980809825369411}