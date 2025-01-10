"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:134937 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a coenzyme A derivative where the acyl group is a fatty acid chain 
    with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for thioester linkage (C(=O)S)
    thioester_smarts = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Check for Coenzyme A moiety by searching for the adenine ring
    adenine_smarts = Chem.MolFromSmarts('n1cnc2c(ncnc12)')  # Adenine ring
    if not mol.HasSubstructMatch(adenine_smarts):
        return False, "Coenzyme A moiety not found"

    # For each thioester linkage, check the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon
        sulfur_idx = match[1]      # Index of sulfur atom

        # Start traversal from the carbonyl carbon
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        visited_idxs = {carbonyl_c_idx, sulfur_idx}

        # Find the alpha carbon (next carbon attached to carbonyl carbon)
        alpha_carbon = None
        for nbr in carbonyl_c.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_idxs:
                alpha_carbon = nbr
                break
        if alpha_carbon is None:
            continue  # No alpha carbon found

        visited_idxs.add(alpha_carbon.GetIdx())

        # Find the beta carbon (next carbon after alpha carbon)
        beta_carbon = None
        for nbr in alpha_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_idxs:
                beta_carbon = nbr
                break
        if beta_carbon is None:
            continue  # No beta carbon found

        visited_idxs.add(beta_carbon.GetIdx())

        # Find the gamma carbon (third carbon after carbonyl carbon)
        gamma_carbon = None
        for nbr in beta_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_idxs:
                gamma_carbon = nbr
                break
        if gamma_carbon is None:
            continue  # No gamma carbon found

        visited_idxs.add(gamma_carbon.GetIdx())

        # Check if gamma carbon has an OH group attached
        has_OH = False
        for nbr in gamma_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen atom
                if nbr.GetDegree() == 1:  # Hydroxyl group (-OH)
                    has_OH = True
                    break
        if not has_OH:
            continue  # No hydroxyl group at position 3

        # Optionally, verify that the acyl chain is sufficiently long (e.g., at least 4 carbons)
        chain_length = 3  # Already have alpha, beta, gamma carbons
        current_atom = gamma_carbon
        while True:
            next_carbon = None
            for nbr in current_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_idxs:
                    next_carbon = nbr
                    break
            if next_carbon is None:
                break  # End of chain
            chain_length += 1
            visited_idxs.add(next_carbon.GetIdx())
            current_atom = next_carbon

        if chain_length < 4:
            continue  # Acyl chain too short to be a fatty acid

        # All conditions met
        return True, "Contains 3-hydroxy fatty acyl-CoA structure"

    # No matching structures found
    return False, "Does not match 3-hydroxy fatty acyl-CoA structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:134937',
        'name': '3-hydroxy fatty acyl-CoA',
        'definition': 'A hydroxy fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.',
        'parents': ['CHEBI:64479', 'CHEBI:83875']
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
    }
}