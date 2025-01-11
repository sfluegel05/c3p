"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (VLCFA) based on its SMILES string.
    VLCFAs are characterized by having a carbon chain longer than 20 carbons with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VLCFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group -O=C(O) pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid found"
        
    # Count carbons in the longest continuous carbon chain
    longest_chain_length = max(len(chain) for chain in rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False) if isinstance(chain, tuple))
    if longest_chain_length <= 20:
        return False, f"Longest carbon chain is {longest_chain_length} carbons, VLCFAs require more than 20"

    return True, "Contains a very long carbon chain with carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27283',
                          'name': 'very long-chain fatty acid',
                          'definition': 'A fatty acid which has a chain length '
                                        'greater than C22. Very long-chain '
                                        'fatty acids which have a chain length '
                                        'greater than C27 are also known as '
                                        'ultra-long-chain fatty acids.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': ['LIPID_MAPS_instance:LMFA01010025'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Longest carbon chain is 15 '
                                               'carbons, VLCFAs require more '
                                               'than 20'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Longest carbon chain is 13 '
                                               'carbons, VLCFAs require more '
                                               'than 20'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': 'C1C[C@@H]([C@@H](O[C@@H]1CCNC(=O)C2CCOCC2)CO)NC(=O)NC3=CC(=CC=C3)Cl',
                                     'name': 'N-[2-[(2S,5S,6R)-5-[[(3-chloroanilino)-oxomethyl]amino]-6-(hydroxymethyl)-2-oxanyl]ethyl]-4-oxanecarboxamide',
                                     'reason': 'No carboxylic acid found'},
                                 {   'smiles': '[H][C@@]1(CNc2nc(N)[nH]c(=O)c2N1)[C@@H](O)[C@H](C)O',
                                     'name': 'sapropterin',
                                     'reason': 'No carboxylic acid found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify very long-chain fatty acids (VLCFAs) using their '
                 "SMILES strings, it's crucial to define the structural "
                 'pattern they exhibit. VLCFAs are characterized by having '
                 'long carbon chains typically over 20 carbons in length, '
                 'often with carboxylic acid end groups. Additionally, they '
                 'may contain unsaturations or functional groups like hydroxy '
                 'groups.\n'
                 '\n'
                 "Here's a structured approach to implementing such a "
                 'classifier:\n'
                 '\n'
                 '1. Parse the SMILES string to validate it as a chemically '
                 'valid structure.\n'
                 '2. Identify the presence of a carboxylic acid group.\n'
                 '3. Count the number of carbon atoms in the longest '
                 'continuous carbon chain.\n'
                 '4. Check for the presence of possible unsaturations '
                 '(double/triple bonds).\n'
                 '5. Ensure the carbon chain length exceeds a certain '
                 'threshold (usually over 20 carbons for a molecule to be '
                 'considered a VLCFA).\n'
                 '6. Return the result with reasoning.\n'
                 '\n'
                 'We should focus on confirming the core characteristic of '
                 'VLCFAs: their exceptionally long hydrocarbon chain with a '
                 'carboxylic acid functional group.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 201,
    'num_false_positives': 100,
    'num_true_negatives': 210,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6677740863787376,
    'recall': 1.0,
    'f1': 0.800796812749004,
    'accuracy': 0.8043052837573386,
    'negative_predictive_value': 1.0}