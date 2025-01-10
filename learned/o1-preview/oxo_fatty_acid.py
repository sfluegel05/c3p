"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:52214 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for aldehyde group (excluding carboxylic acid)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    # Check for ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Exclude aldehyde and ketone groups that are part of the carboxylic acid
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_acid_atoms = set()
    for match in carboxylic_acid_matches:
        carboxylic_acid_atoms.update(match)

    valid_oxo_groups = []
    for match in aldehyde_matches + ketone_matches:
        # Check if any atoms in the oxo group overlap with carboxylic acid atoms
        if not any(atom_idx in carboxylic_acid_atoms for atom_idx in match):
            valid_oxo_groups.append(match)

    if not valid_oxo_groups:
        return False, "No aldehyde or ketone group found in addition to carboxylic acid"

    # Check for long aliphatic chain (fatty acid chain)
    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 8:
        return False, f"Too few carbon atoms ({num_carbons}) for a fatty acid"

    # Check for presence of long carbon chain
    # Find the longest carbon chain in the molecule
    chains = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    max_chain_length = 0
    for chain in chains:
        if all(atom.GetAtomicNum() in [6,1] for atom in chain.GetAtoms()):
            length = chain.GetNumHeavyAtoms()
            if length > max_chain_length:
                max_chain_length = length

    if max_chain_length < 8:
        return False, f"No long aliphatic chain found (maximum chain length: {max_chain_length})"

    return True, "Contains carboxylic acid and aldehyde/ketone groups with a long aliphatic chain"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52214',
        'name': 'oxo fatty acid',
        'definition': 'Any fatty acid containing at least one aldehydic or ketonic group in addition to the carboxylic acid group.',
        'parents': ['CHEBI:25360', 'CHEBI:25352']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}