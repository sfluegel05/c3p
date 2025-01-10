"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as 'Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid carbon atom
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    matches = mol.GetSubstructMatches(carboxy_pattern)
    if not matches:
        return False, "No carboxylic acid group found"

    carboxylic_carbon_idx = matches[0][0]  # Carbon atom of the carboxylic group

    # Traverse the carbon chain starting from the carboxylic carbon
    visited = set()
    stack = [(carboxylic_carbon_idx, -1)]  # (current_atom_idx, parent_atom_idx)

    chain_atom_indices = []
    branching = False

    while stack:
        current_idx, parent_idx = stack.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Only consider carbon atoms for the main chain
        chain_atom_indices.append(current_idx)

        # Get neighboring carbon atoms excluding the parent atom
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != parent_idx]
        num_carbon_neighbors = len(neighbors)

        if num_carbon_neighbors > 1:
            branching = True  # More than one carbon neighbor (excluding parent) indicates branching

        for nbr in neighbors:
            nbr_idx = nbr.GetIdx()
            stack.append((nbr_idx, current_idx))

    if branching:
        return False, "Branching detected in the carbon chain"

    carbon_count = len(chain_atom_indices)
    if carbon_count != 18:
        return False, f"Main carbon chain has {carbon_count} carbons, expected 18"

    # Count C=C double bonds between carbons in the main chain
    double_bond_count = 0
    bond_indices = set()
    for i in range(len(chain_atom_indices) - 1):
        atom_idx1 = chain_atom_indices[i]
        atom_idx2 = chain_atom_indices[i + 1]
        bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is not None:
            bond_indices.add(bond.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1

    if double_bond_count != 2:
        return False, f"Found {double_bond_count} C=C double bonds in the main chain, expected 2"

    return True, "Molecule is a straight-chain C18 fatty acid with two C=C double bonds"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'octadecadienoic acid',
        'definition': 'Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.',
        'parents': None
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