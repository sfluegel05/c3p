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

    # Identify carboxylic acid group
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if not carboxy_matches:
        return False, "No terminal carboxylic acid group found"

    # Get the carboxylic carbon atom index
    carboxylic_carbon_idx = carboxy_matches[0][0]

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

        # Exclude atoms that are not carbon
        if atom.GetAtomicNum() != 6:
            # Allow oxygen atoms only if they are part of the carboxylic group
            if current_idx in [idx for match in carboxy_matches for idx in match]:
                continue
            else:
                return False, f"Atom {current_idx} is not carbon (found {atom.GetSymbol()}) in the main chain"

        chain_atom_indices.append(current_idx)

        # Get neighboring carbon atoms excluding the parent atom
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != parent_idx]

        # Check for branching (more than 2 neighbors means branching in a linear chain)
        carbon_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            branching = True

        for nbr in neighbors:
            nbr_idx = nbr.GetIdx()
            nbr_atom = mol.GetAtomWithIdx(nbr_idx)
            # Exclude non-carbon atoms (heteroatoms) in the chain
            if nbr_atom.GetAtomicNum() != 6:
                # Allow oxygen atoms only if they are part of the carboxylic group
                if nbr_idx in [idx for match in carboxy_matches for idx in match]:
                    continue
                else:
                    return False, f"Non-carbon atom {nbr_atom.GetSymbol()} found in the main chain"
            if nbr_idx != parent_idx:
                stack.append((nbr_idx, current_idx))

    if branching:
        return False, "Branching detected in the carbon chain"

    carbon_count = len(chain_atom_indices)
    if carbon_count != 18:
        return False, f"Main carbon chain has {carbon_count} carbons, expected 18"

    # Check for triple bonds or rings
    ring_info = mol.GetRingInfo()
    if ring_info.IsBondInRingOfSize(3):
        return False, "Cyclopropane ring detected"
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            return False, "Triple bond detected in the main chain"

    # Count C=C double bonds between carbons in the main chain
    double_bond_count = 0
    for i in range(len(chain_atom_indices) - 1):
        atom_idx1 = chain_atom_indices[i]
        atom_idx2 = chain_atom_indices[i + 1]
        bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is not None:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_count += 1
            elif bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                return False, "Non-single/double bond detected in the main chain"

    if double_bond_count != 2:
        return False, f"Found {double_bond_count} C=C double bonds in the main chain, expected 2"

    # Check for additional functional groups (other than the carboxyl group)
    heteroatom_pattern = Chem.MolFromSmarts("[!C;!H]")
    heteroatom_matches = mol.GetSubstructMatches(heteroatom_pattern)
    allowed_heteroatoms = set()
    for match in carboxy_matches:
        allowed_heteroatoms.update(match)
    for idx in [idx for match in heteroatom_matches for idx in match]:
        if idx not in allowed_heteroatoms:
            atom = mol.GetAtomWithIdx(idx)
            return False, f"Additional functional group detected: {atom.GetSymbol()} at atom index {idx}"

    return True, "Molecule is a straight-chain C18 fatty acid with two C=C double bonds and no additional functional groups"

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
    'attempt': 3,
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