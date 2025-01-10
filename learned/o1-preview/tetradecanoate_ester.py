"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester formed from tetradecanoic acid (myristic acid) and an alcohol.
    The acyl chain must contain exactly 14 carbons (including the carbonyl carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester SMARTS pattern
    ester_pattern = Chem.MolFromSmarts("[C](=O)[O][#6]")  # Ester group with carbonyl carbon connected to oxygen and carbon
    if ester_pattern is None:
        return False, "Failed to create ester SMARTS pattern"

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if acyl chain is tetradecanoyl (14 carbons)
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon
        ester_o_idx = match[2]     # Index of ester oxygen

        # Initialize variables for traversal
        visited = set()
        chain_length = 0
        is_linear = True
        is_saturated = True

        # Traverse acyl chain starting from carbonyl carbon (excluding ester oxygen side)
        def traverse_acyl_chain(atom_idx, prev_atom_idx):
            nonlocal chain_length, is_linear, is_saturated

            atom = mol.GetAtomWithIdx(atom_idx)
            if atom_idx in visited:
                return
            visited.add(atom_idx)

            if atom.GetAtomicNum() != 6:
                return  # Non-carbon atom

            # Increment chain length
            chain_length += 1

            # Check for branching (more than two carbon neighbors, excluding previous atom)
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom_idx]
            if len(neighbor_carbons) > 1:
                is_linear = False

            # Check for unsaturation (non-single bonds)
            for bond in atom.GetBonds():
                if bond.GetBeginAtomIdx() == prev_atom_idx or bond.GetEndAtomIdx() == prev_atom_idx:
                    continue  # Skip bond to previous atom
                if bond.GetBondTypeAsDouble() != 1.0:
                    is_saturated = False

            # Continue traversal to neighbor carbons (excluding previous atom and ester oxygen)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == prev_atom_idx or neighbor_idx == ester_o_idx:
                    continue
                if neighbor.GetAtomicNum() == 6:
                    traverse_acyl_chain(neighbor_idx, atom_idx)

        # Start traversal from carbonyl carbon
        traverse_acyl_chain(carbonyl_c_idx, ester_o_idx)

        # Check if acyl chain has exactly 14 carbons
        if chain_length == 14 and is_linear and is_saturated:
            return True, "Contains tetradecanoate ester group"

    return False, "Does not contain tetradecanoate ester group"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'tetradecanoate ester',
        'definition': 'A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid (also known as myristic acid) with a hydroxy group of an alcohol or phenol.',
        'parents': []
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