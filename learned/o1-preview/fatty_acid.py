"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:28868 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons,
    which may be saturated or unsaturated, and can include hydroxyl groups or minor substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule is not aromatic
    if mol.HasSubstructMatch(Chem.MolFromSmarts('a')):
        return False, "Molecule contains aromatic rings, which is not characteristic of fatty acids"

    # Identify carboxylic acid groups: [C](=O)[O][H]
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acids = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acids) != 1:
        return False, f"Found {len(carboxylic_acids)} carboxylic acid groups, need exactly 1"

    # Get the carboxyl carbon atom index
    carboxyl_c_index = carboxylic_acids[0][0]

    # Allowed atoms in fatty acids
    allowed_atomic_nums = {1, 6, 7, 8, 9, 16, 17, 35, 53}  # H, C, N, O, F, S, Cl, Br, I

    # Check for disallowed elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"

    # Traverse the molecule to find the longest carbon chain connected to the carboxyl carbon
    visited = set()

    def traverse_chain(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 6:  # Only consider carbon atoms
            return []
        visited.add(atom_idx)
        chain = [atom_idx]
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:
                chain.extend(traverse_chain(n_idx))
        return chain

    # Start traversal from the carbon next to the carboxyl carbon
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_c_index)
    chains = []
    for neighbor in carboxyl_carbon.GetNeighbors():
        n_idx = neighbor.GetIdx()
        if mol.GetAtomWithIdx(n_idx).GetAtomicNum() == 6:
            visited.clear()
            chain = traverse_chain(n_idx)
            chains.append(chain)

    if not chains:
        return False, "No carbon chain found attached to the carboxyl group"

    # Find the longest carbon chain
    longest_chain = max(chains, key=len)
    chain_length = len(set(longest_chain))  # Number of unique carbons in the chain

    # Check if chain length is within typical fatty acid range
    if chain_length < 4:
        return False, f"Chain length {chain_length} is too short for a fatty acid"
    elif chain_length > 28:
        return False, f"Chain length {chain_length} is too long for typical fatty acids"

    # Ensure the chain is aliphatic (no rings)
    for idx in longest_chain:
        atom = mol.GetAtomWithIdx(idx)
        if atom.IsInRing():
            return False, "Carbon chain contains ring structures, which is not characteristic of fatty acids"

    return True, f"Contains a single carboxylic acid group and an aliphatic chain of {chain_length} carbons"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:28868',
        'name': 'fatty acid',
        'definition': 'Any aliphatic monocarboxylic acid derived from or contained in esterified form in an animal or vegetable fat, oil or wax. Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched and even-numbered), which may be saturated or unsaturated. By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.',
        'parents': ['CHEBI:25384', 'CHEBI:17436']
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
    'attempt': 1,
    'success': False,
    'best': False,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0
}