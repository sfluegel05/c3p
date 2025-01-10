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

    # Check that the molecule has exactly one carboxylic acid group (C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxy_matches) != 1:
        return False, "Molecule does not have exactly one carboxylic acid group"
    
    # Check that the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures"

    # Check that the molecule is aliphatic (no aromatic atoms)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms"

    # Check that the molecule contains only C, H, O atoms
    allowed_atoms = {6, 1, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Molecule contains forbidden atom type: {atom.GetSymbol()}"

    # Get the carboxyl carbon atom index
    carboxyl_carbon_idx = carboxy_matches[0][0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Get the carbon connected to the carboxyl carbon (the alpha carbon)
    neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(neighbors) != 1:
        return False, "Carboxyl carbon does not have exactly one carbon neighbor"
    alpha_carbon = neighbors[0]
    alpha_carbon_idx = alpha_carbon.GetIdx()

    # Traverse the hydrocarbon chain starting from the alpha carbon
    visited = set()
    branching_found = False
    chain_length = 0

    from collections import deque
    queue = deque()
    queue.append((alpha_carbon_idx, None))  # (current atom idx, previous atom idx)

    while queue:
        current_idx, parent_idx = queue.popleft()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            return False, "Hydrocarbon chain contains non-carbon atoms"
        chain_length += 1
        neighbor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() != parent_idx]
        carbon_neighbors = [idx for idx in neighbor_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        # Check for branching (more than one carbon neighbor besides parent)
        if len(carbon_neighbors) > 1:
            branching_found = True
        # Add neighbors to the queue
        for nbr_idx in carbon_neighbors:
            queue.append((nbr_idx, current_idx))

    # Ensure the chain is of sufficient length (e.g., at least 8 carbons)
    if chain_length < 8:
        return False, f"Hydrocarbon chain too short ({chain_length} carbons), not a fatty acid"

    if branching_found:
        return True, "Branched-chain fatty acid detected"

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}