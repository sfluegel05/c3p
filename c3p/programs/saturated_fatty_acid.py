"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: saturated fatty acids
Definition: Any fatty acid containing no carbon to carbon multiple bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Check for absence of carbon-carbon multiple bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Contains carbon-carbon double bonds"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds"
    
    # Count carbons and check chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, "Carbon chain too short for fatty acid"
        
    # Check atom types (allowing for deuterium)
    allowed_atoms = {1, 6, 8}  # H, C, O
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains atoms other than C, H, O"
        
    # Check that carboxylic acid is terminal
    # First get the carbon of the carboxyl group
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_matches[0][0])
    
    # Count non-oxygen connections to carboxyl carbon
    non_oxygen_connections = sum(1 for neighbor in carboxyl_carbon.GetNeighbors() 
                               if neighbor.GetAtomicNum() != 8)
    
    if non_oxygen_connections > 1:
        return False, "Carboxylic acid group is not terminal"
        
    # Count oxygens - should have exactly 2 for the carboxyl group
    # (plus any additional OH groups in hydroxy fatty acids)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for carboxylic acid"
        
    return True, "Saturated fatty acid with terminal carboxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26607',
                          'name': 'saturated fatty acid',
                          'definition': 'Any fatty acid containing no carbon '
                                        'to carbon multiple bonds. Known to '
                                        'produce adverse biological effects '
                                        'when ingested to excess.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': [   'PMID:16492686',
                                       'PMID:19763019',
                                       'PMID:20237329'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Multiple carboxylic acid '
                                               'groups found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Contains carbon-carbon double '
                                               'bonds'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC)COC(CCCCCCCCCCCCCCC)=O',
                                      'name': '1-hexadecanoyl-2-methyl-sn-glycero-3-phosphocholine',
                                      'reason': 'No carboxylic acid group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify saturated fatty '
                 "acids. Here's my thought process:\n"
                 '\n'
                 '1. First, we need to check if the molecule is a carboxylic '
                 'acid by looking for the -COOH group\n'
                 '2. Then we need to verify it has no double or triple bonds '
                 'between carbons\n'
                 '3. We should check that it has a reasonable carbon chain '
                 'length (typically >3 carbons)\n'
                 '4. The molecule should be composed of only C, H, and O atoms '
                 '(with possible exceptions for isotopes like deuterium)\n'
                 '5. The carboxylic acid should be terminal (at the end of the '
                 'chain)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 89,
    'num_false_positives': 100,
    'num_true_negatives': 5002,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.4708994708994709,
    'recall': 0.9888888888888889,
    'f1': 0.6379928315412186,
    'accuracy': 0.9805469953775039,
    'negative_predictive_value': 0.9998001199280432}