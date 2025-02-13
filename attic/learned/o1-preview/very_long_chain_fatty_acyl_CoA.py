"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:51953 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester group pattern: C(=O)S
    thioester_smarts = 'C(=O)S'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Invalid thioester SMARTS pattern"

    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Simplified CoA moiety pattern
    coa_smarts = 'NC(=O)CCNC(=O)C(O)C(C)(C)CO'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # For each thioester linkage, attempt to find the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon
        sulfur_idx = match[2]      # Index of sulfur atom

        # Get the atom connected to carbonyl carbon that is not sulfur
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_c_idx)
        acyl_chain_atom = None
        for neighbor in carbonyl_carbon.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != sulfur_idx:
                acyl_chain_atom = neighbor
                break

        if acyl_chain_atom is None:
            continue  # Move to next thioester linkage if no acyl chain found

        # Traverse the acyl chain starting from acyl_chain_atom
        acyl_chain_carbons = set()
        atoms_to_visit = [(acyl_chain_atom.GetIdx(), carbonyl_c_idx)]  # (current_atom_idx, previous_atom_idx)
        while atoms_to_visit:
            current_atom_idx, previous_atom_idx = atoms_to_visit.pop()
            current_atom = mol.GetAtomWithIdx(current_atom_idx)
            if current_atom.GetAtomicNum() != 6:
                continue  # Only consider carbon atoms
            acyl_chain_carbons.add(current_atom_idx)
            for neighbor in current_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == previous_atom_idx:
                    continue  # Avoid going back to previous atom
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neighbor_atom.GetAtomicNum() == 6 and neighbor_idx not in acyl_chain_carbons:
                    atoms_to_visit.append((neighbor_idx, current_atom_idx))

        chain_length = len(acyl_chain_carbons)
        if chain_length > 22:
            return True, f"Acyl chain length is {chain_length}, which is greater than 22 carbons"
        else:
            return False, f"Acyl chain length is {chain_length}, which is not greater than 22 carbons"

    return False, "No valid acyl chain found attached to CoA"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:51953',
        'name': 'very long-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.',
        'parents': ['CHEBI:37554']},
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
        'test_proportion': 0.1},
    'message': None,
    'attempt': 3
}