"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    A branched-chain fatty acid is a fatty acid in which the hydrocarbon chain has one or more
    alkyl substituents (branches). The fatty acyl chain is usually saturated and the substituent
    is a methyl group; however, unsaturated BCFAs are found in marine animals, and branches
    other than methyl are found in microbial lipids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Get the carboxyl carbon atom index
    carboxyl_carbon_idx = carboxy_matches[0][0]

    # Find all acyclic paths (i.e., chains) starting from the carboxyl carbon
    # We will identify the longest carbon chain as the main fatty acyl chain
    from collections import deque

    visited = set()
    max_chain = []
    branching_found = False

    queue = deque()
    queue.append((carboxyl_carbon_idx, [carboxyl_carbon_idx]))

    while queue:
        current_atom_idx, path = queue.popleft()
        visited.add(current_atom_idx)
        current_atom = mol.GetAtomWithIdx(current_atom_idx)

        # Track the longest chain
        if len(path) > len(max_chain):
            max_chain = path

        # Get neighbors excluding oxygens (to stay within the carbon chain)
        neighbors = [
            nbr.GetIdx() for nbr in current_atom.GetNeighbors()
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in path
        ]

        # Check for branching: if current atom has more than one carbon neighbor not in path
        num_carbon_neighbors = sum(1 for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
        if num_carbon_neighbors > 2:
            branching_found = True

        for nbr_idx in neighbors:
            queue.append((nbr_idx, path + [nbr_idx]))

    if not max_chain:
        return False, "No hydrocarbon chain found"

    # After finding the main chain, check for branches
    main_chain_set = set(max_chain)
    for atom_idx in max_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6 and nbr_idx not in main_chain_set:
                # Found a carbon branch off the main chain
                return True, "Branched chain found in hydrocarbon chain"

    if branching_found:
        return True, "Branched chain found in hydrocarbon chain"

    return False, "No branching in hydrocarbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25941',
        'name': 'branched-chain fatty acid',
        'definition': 'Any fatty acid in which the parent hydrocarbon chain has one or more alkyl substituents; a common component in animal and bacterial lipids. The fatty acyl chain is usually saturated and the substituent a methyl group; however, unsaturated BCFAs are found in marine animals, and branches other than methyl are found in microbial lipids.',
        'parents': ['CHEBI:15740', 'CHEBI:25359']
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
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}