"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:24863 organoiodine compound
Definition: A compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound must contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    
    if not iodine_atoms:
        return False, "No iodine atoms found"

    # Check if any iodine is bonded to carbon
    for iodine in iodine_atoms:
        for neighbor in iodine.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                return True, "Contains at least one carbon-iodine bond"
    
    return False, "No carbon-iodine bonds found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37142',
                          'name': 'organoiodine compound',
                          'definition': 'An organoiodine compound is a '
                                        'compound containing at least one '
                                        'carbon-iodine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:24860'],
                          'xrefs': [   'MetaCyc:Organoiodine-Compounds',
                                       'Wikipedia:Organoiodine_compound'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No iodine atoms found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No iodine atoms found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': "I'll create a program to classify organoiodine compounds. "
                 'The key requirement is to detect at least one carbon-iodine '
                 "bond in the molecule. Here's my solution:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 99,
    'num_false_positives': 100,
    'num_true_negatives': 73352,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.49748743718592964,
    'recall': 1.0,
    'f1': 0.6644295302013422,
    'accuracy': 0.9986403991788011,
    'negative_predictive_value': 1.0}