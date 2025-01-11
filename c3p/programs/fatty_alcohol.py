"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain of 3 to >27 atoms, which may be saturated/unsaturated and branched/unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for at least one hydroxyl (-OH) group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group pattern
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check carbon count range
    if carbon_count < 3:
        return False, f"Too few carbon atoms ({carbon_count}), need at least 3"
    # While unspecified max, more than 27 is in the description, no check needed unless an upper threshold is found relevant

    return True, "Contains a carbon chain with a hydroxyl group(s), characteristic of fatty alcohols"

# The function is designed to identify key features that generally define fatty alcohols as per the defined chemical class.


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24026',
                          'name': 'fatty alcohol',
                          'definition': 'An aliphatic alcohol consisting of a '
                                        'chain of 3 to greater than 27 carbon '
                                        'atoms. Fatty alcohols may be '
                                        'saturated or unsaturated and may be '
                                        'branched or unbranched.',
                          'parents': ['CHEBI:30879', 'CHEBI:61697'],
                          'xrefs': [   'LIPID_MAPS_class:LMFA05',
                                       'MetaCyc:Fatty-Alcohols',
                                       'Wikipedia:Fatty_alcohol'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'C12=C(CCC[C@@]1(CC[C@H](C2)C(=C)CCC=C(C)C)C)C',
                                     'name': 'Clavulatriene A',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'CCOC(=O)n1ccn(C)c1=S',
                                     'name': 'carbimazole',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'C1CC(=NC1)NNC(=O)CNC2=CC3=CC=CC=C3C=C2',
                                     'name': "N'-(3,4-dihydro-2H-pyrrol-5-yl)-2-(2-naphthalenylamino)acetohydrazide",
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'NC(=N)c1ccc(cc1)-c1cc2ccc(cc2[nH]1)C(N)=N',
                                     'name': 'DAPI',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'C1CC(CCC1S(=O)(=O)C2=CC3=C(C=C2)OCCO3)S(=O)(=O)C4=C(C=CC=C4F)F',
                                     'name': '6-[4-(2,6-difluorophenyl)sulfonylcyclohexyl]sulfonyl-2,3-dihydro-1,4-benzodioxin',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'CN(C(=O)C1=CC2=C(OC(F)(F)O2)C=C1)C1=CC=CC(C(=O)NC2=C(Br)C=C(C=C2OC(F)F)C(F)(C(F)(F)F)C(F)(F)F)=C1F',
                                     'name': 'piperflanilide',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'COC1(CCN(CC1)CCC2=CNC3=C2C=C(C=C3)F)C[S@@](=O)C4=CC=CC=C4',
                                     'name': '5-fluoro-3-[2-[4-methoxy-4-[[(R)-phenylsulfinyl]methyl]-1-piperidinyl]ethyl]-1H-indole',
                                     'reason': 'No hydroxyl group found'},
                                 {   'smiles': 'C1=CC=C2C(=C1)C=3CCN4C(=C5C=CC=6C=NC=C(C65)C4=O)C3N2',
                                     'name': 'Naulafine',
                                     'reason': 'No hydroxyl group found'}],
    'sample_false_negatives': [   {   'smiles': 'Clc1c2COC(=O)c2c(Cl)c(Cl)c1Cl',
                                      'name': '4,5,6,7-tetrachloro-2-benzofuran-1(3H)-one',
                                      'reason': 'No hydroxyl group found'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical as a fatty alcohol based on its '
                 'SMILES string, we need to implement a function that checks '
                 'the structural features of fatty alcohols:\n'
                 '\n'
                 '1. **Carbon Atom Range**: Ensure the molecule contains a '
                 'carbon chain of at least 3 carbon atoms and can extend to '
                 'more than 27 carbon atoms.\n'
                 '\n'
                 '2. **Alcohol Functional Group**: The presence of one or more '
                 'hydroxyl groups (-OH) attached to the carbon chain.\n'
                 '\n'
                 '3. **Saturation and Branching**: The chain can be saturated '
                 '(single bonds) or unsaturated (double or triple bonds), and '
                 'it can be branched or unbranched. However, the primary '
                 'feature remains the presence of a chain of hydrocarbons and '
                 'hydroxyl groups.\n'
                 '\n'
                 'We will use RDKit to parse the SMILES string and identify '
                 'these structural components.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 433,
    'num_false_positives': 100,
    'num_true_negatives': 51,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.8123827392120075,
    'recall': 0.9976958525345622,
    'f1': 0.8955532574974147,
    'accuracy': 0.8273504273504273,
    'negative_predictive_value': 0.9807692307692307}