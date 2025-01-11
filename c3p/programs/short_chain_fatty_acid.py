"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:27283 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of carboxylic acid groups - should be exactly one
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches > 1:
        return False, "Multiple carboxylic acid groups found"
        
    # Check for aromatic rings - should have none
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
        
    # Count carbons (excluding the carboxyl carbon)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    chain_carbons = carbon_count - 1  # Subtract the carboxyl carbon
    
    if chain_carbons >= 6:
        return False, f"Chain too long ({chain_carbons} carbons excluding carboxyl group)"
    if chain_carbons < 1:
        return False, "No carbon chain found"
        
    # Check for non-hydrocarbon atoms (except the carboxylic acid oxygens)
    allowed_atoms = {6, 1, 8}  # C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains non-hydrocarbon atoms other than carboxylic acid"
            
    # Count oxygens - should be exactly 2 (from carboxylic acid)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 2:
        return False, "Contains oxygen atoms other than carboxylic acid"
        
    # All checks passed
    return True, f"Aliphatic monocarboxylic acid with {chain_carbons} carbons in chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27283',
        'name': 'short-chain fatty acid',
        'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. '
                     'If any non-hydrocarbon substituent is present, the compound is not normally '
                     'regarded as a short-chain fatty acid.',
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26666',
                          'name': 'short-chain fatty acid',
                          'definition': 'An aliphatic monocarboxylic acid with '
                                        'a chain length of less than C6. If '
                                        'any non-hydrocarbon substituent is '
                                        'present, the compound is not normally '
                                        'regarded as a short-chain fatty acid.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': [   'PMID:16633129',
                                       'PMID:16870803',
                                       'PMID:18203540',
                                       'PMID:20148677'],
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
                                     'reason': 'Chain too long (21 carbons '
                                               'excluding carboxyl group)'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'OC(C(CC)C)C(O)=O',
                                      'name': '2-hydroxy-3-methylpentanoic '
                                              'acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'CCC(=O)C(O)=O',
                                      'name': '2-oxobutanoic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'OC(C[C@H](CC)O)=O',
                                      'name': '(S)-3-hydroxypentanoic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'CCC(CO)C(O)=O',
                                      'name': '2-ethylhydracrylic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': '[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C',
                                      'name': 'heliosupine',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(CC(CC)O)=O',
                                      'name': '3-hydroxypentanoic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'OC(C[C@@H](CC)O)=O',
                                      'name': '(R)-3-hydroxypentanoic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'OCCC(O)=O',
                                      'name': '3-hydroxypropionic acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'},
                                  {   'smiles': 'CCC/C(=C/C(O)=O)/C',
                                      'name': '(2E)-3-methylhex-2-enoic acid',
                                      'reason': 'Chain too long (6 carbons '
                                                'excluding carboxyl group)'},
                                  {   'smiles': 'O[C@H]([C@@H](CC)C)C(O)=O',
                                      'name': '(2R,3R)-2-hydroxy-3-methylpentanoic '
                                              'acid',
                                      'reason': 'Contains oxygen atoms other '
                                                'than carboxylic acid'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify short-chain fatty '
                 'acids. Based on the definition and examples, we need to '
                 'check for:\n'
                 '\n'
                 '1. Presence of a carboxylic acid group\n'
                 '2. Chain length less than C6\n'
                 '3. No non-hydrocarbon substituents (except the carboxylic '
                 'acid group)\n'
                 '4. Aliphatic structure (no aromatic rings)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 29,
    'num_false_positives': 31,
    'num_true_negatives': 142230,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.48333333333333334,
    'recall': 0.7435897435897436,
    'f1': 0.5858585858585859,
    'accuracy': 0.9997118763176388,
    'negative_predictive_value': 0.999929696287964}