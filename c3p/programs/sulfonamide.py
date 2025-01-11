"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide contains the functional group RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sulfonamide pattern
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)[N]")
    
    # Check for the sulfonamide group
    if not mol.HasSubstructMatch(sulfonamide_pattern):
        return False, "No sulfonamide group found"
    
    return True, "Contains a sulfonamide group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35358',
                          'name': 'sulfonamide',
                          'definition': 'An amide of a sulfonic acid '
                                        "RS(=O)2NR'2.",
                          'parents': [   'CHEBI:33256',
                                         'CHEBI:33261',
                                         'CHEBI:33552'],
                          'xrefs': [   'PMID:11498380',
                                       'PMID:2434548',
                                       'PMID:26811268',
                                       'PMID:26832216',
                                       'Wikipedia:Sulfonamide'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No sulfonamide group found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No sulfonamide group found'}],
    'sample_false_negatives': [   {   'smiles': 'COC1=CC(=C(C=C1Cl)OC)S(=O)(=O)N2C3=CC=CC=C3N=N2',
                                      'name': '1-(4-chloro-2,5-dimethoxyphenyl)sulfonylbenzotriazole',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'CCOC1=C(C=CC(=C1)S(=O)(=O)N2C(=CC(=N2)C)C)Br',
                                      'name': '1-(4-bromo-3-ethoxyphenyl)sulfonyl-3,5-dimethylpyrazole',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'C1=CC=C(C=C1)C2=NC=C(S2)C3=NN(C=C3)S(=O)(=O)C4=CC=CC=C4',
                                      'name': '5-[1-(benzenesulfonyl)-3-pyrazolyl]-2-phenylthiazole',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'CC1=CC=C(C=C1)S(=O)(=O)N2C3=CC=CC=C3N=C2N',
                                      'name': '1-(4-methylphenyl)sulfonyl-2-benzimidazolamine',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'CC1=NN(C(=C1S(=O)(=O)C2=CC=CC=C2)OC(=O)C3CCCC3)C(C)(C)C',
                                      'name': 'cyclopentanecarboxylic acid '
                                              '[4-(benzenesulfonyl)-2-tert-butyl-5-methyl-3-pyrazolyl] '
                                              'ester',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'C1CCN(CC1)CC(CS(=O)(=O)C2=CC=CC=C2)N3CCCCC3',
                                      'name': '1-[3-(benzenesulfonyl)-2-(1-piperidinyl)propyl]piperidine',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'C1C(=O)NC(=CS(=O)(=O)C2=CC=C(C=C2)Cl)S1',
                                      'name': '2-[(4-chlorophenyl)sulfonylmethylidene]-4-thiazolidinone',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'C1=CC(=CC=C1S(C=2C=CC(=CC2)NC(S(O)(=O)=O)[C@@H]([C@@H](O)[C@@H]([C@H](O)CO)O)O)(=O)=O)NC([C@@H]([C@@H](O)[C@@H]([C@H](O)CO)O)O)S(O)(=O)=O',
                                      'name': 'glucosulfone',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'CS(=O)(=O)C1=CC=C(C=C1)C([C@H](CF)NC(=O)C(Cl)Cl)O',
                                      'name': '2,2-dichloro-N-[(2R)-3-fluoro-1-hydroxy-1-(4-methylsulfonylphenyl)propan-2-yl]acetamide',
                                      'reason': 'No sulfonamide group found'},
                                  {   'smiles': 'C1=CC=C2C(=C1)N=C(N2S(=O)(=O)C3=CC=C(C=C3)Cl)N',
                                      'name': '1-(4-chlorophenyl)sulfonyl-2-benzimidazolamine',
                                      'reason': 'No sulfonamide group found'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a sulfonamide, we will check '
                 'for the presence of the functional group sulfonamide, which '
                 'can be identified by the presence of the sulfonyl functional '
                 'group (SO2) attached to a nitrogen atom. The SMILES pattern '
                 'to recognize a sulfonamide is `S(=O)(=O)[N]`.\n'
                 '\n'
                 'Here is a step-by-step breakdown of the approach:\n'
                 '\n'
                 '1. **Import RDKit Libraries**: We need RDKit to parse and '
                 'analyze the SMILES strings.\n'
                 '\n'
                 '2. **Define the Function `is_sulfonamide`**:\n'
                 '   - Parse the input SMILES string into a molecular object '
                 'using RDKit.\n'
                 '   - Define a SMARTS pattern to identify the sulfonamide '
                 'functional group.\n'
                 '   - Check if the molecule contains the specified '
                 'sulfonamide pattern using substructure search.\n'
                 '\n'
                 '3. **Return Results**:\n'
                 '   - If the SMILES contains the sulfonamide pattern, return '
                 '`True` with an appropriate reason.\n'
                 "   - If the SMILES is invalid or doesn't contain the "
                 'sulfonamide pattern, return `False` with a reason explaining '
                 'why it is not a sulfonamide.\n'
                 '\n'
                 "Let's write the code for this:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2290,
    'num_false_positives': 100,
    'num_true_negatives': 2667,
    'num_false_negatives': 40,
    'num_negatives': None,
    'precision': 0.9581589958158996,
    'recall': 0.9828326180257511,
    'f1': 0.9703389830508474,
    'accuracy': 0.9725328624681185,
    'negative_predictive_value': 0.9852234946435168}