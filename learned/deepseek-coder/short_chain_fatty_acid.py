"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a chain length of less than C6.
    If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.

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

    # Check for the presence of a carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check that the molecule is aliphatic (no aromatic rings)
    if mol.GetNumAromaticRings() > 0:
        return False, "Molecule contains aromatic rings"

    # Check for non-hydrocarbon substituents (e.g., halogens, nitrogen, sulfur, etc.)
    non_hydrocarbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]]
    if len(non_hydrocarbon_atoms) > 0:
        return False, "Molecule contains non-hydrocarbon substituents"

    # Calculate the number of carbon atoms in the longest chain
    longest_chain = rdMolDescriptors.CalcLongestChain(mol)
    if longest_chain >= 6:
        return False, f"Chain length is {longest_chain}, which is too long for a short-chain fatty acid"

    # Check that the molecule is a monocarboxylic acid (only one carboxylic acid group)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    return True, "Aliphatic monocarboxylic acid with a chain length of less than C6 and no non-hydrocarbon substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26666',
                          'name': 'short-chain fatty acid',
                          'definition': 'An aliphatic monocarboxylic acid with a chain length of less than C6. If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.',
                          'parents': ['CHEBI:26666', 'CHEBI:26667']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}