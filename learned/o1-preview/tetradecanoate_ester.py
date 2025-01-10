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

    # Define a SMARTS pattern for ester group: carbonyl carbon connected to oxygen connected to carbon
    ester_pattern = Chem.MolFromSmarts("[#6X3](=O)[OX2H0][#6]")
    if ester_pattern is None:
        return False, "Failed to create ester SMARTS pattern"

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain has exactly 14 carbons
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        ester_o_idx = match[1]     # Ester oxygen index
        alcohol_c_idx = match[2]   # Alcohol carbon index

        # Initialize a set to keep track of visited atoms
        visited = set()
        acyl_carbons = set()

        # Traverse the acyl chain starting from the carbonyl carbon (excluding ester oxygen side)
        to_visit = [carbonyl_c_idx]
        while to_visit:
            atom_idx = to_visit.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            if atom.GetAtomicNum() == 6:
                # Add carbon atom to acyl carbons set
                acyl_carbons.add(atom_idx)

                # Visit neighboring atoms
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    # Do not go back to ester oxygen or previously visited atoms
                    if neighbor_idx != ester_o_idx and neighbor_idx not in visited:
                        to_visit.append(neighbor_idx)
            else:
                # Continue traversal for non-carbon atoms (e.g., to handle branching)
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx != ester_o_idx and neighbor_idx not in visited:
                        to_visit.append(neighbor_idx)

        # Count the number of carbons in the acyl chain
        num_carbons = len(acyl_carbons)

        if num_carbons == 14:
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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}