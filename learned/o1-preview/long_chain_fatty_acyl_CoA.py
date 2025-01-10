"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57395 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group of coenzyme A with
    the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thioester linkage pattern S-C(=O)-C
    thioester_smarts = "[#16]-C(=O)-[C]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Failed to construct thioester pattern"

    # Define adenine pattern (part of CoA)
    adenine_smarts = "n1cnc2c1ncnc2N"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Failed to construct adenine pattern"

    # Check if adenine moiety is present (indicates CoA)
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found (CoA moiety missing)"

    # Find thioester linkage(s)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage
    for match in thioester_matches:
        sulfur_idx = match[0]
        carbonyl_c_idx = match[1]
        alpha_c_idx = match[2]

        # Get CoA atom indices by traversing from sulfur (excluding carbonyl carbon)
        CoA_atom_indices = set()
        visited_coa = set()
        def traverse_CoA(atom_idx):
            if atom_idx in visited_coa:
                return
            visited_coa.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx != carbonyl_c_idx and neighbor_idx not in visited_coa:
                    traverse_CoA(neighbor_idx)
        traverse_CoA(sulfur_idx)
        CoA_atom_indices = visited_coa

        # Traverse the fatty acyl chain starting from alpha carbon
        visited_acyl = set()
        carbon_count = 0
        def traverse_acyl_chain(atom_idx):
            nonlocal carbon_count
            if atom_idx in visited_acyl or atom_idx in CoA_atom_indices:
                return
            visited_acyl.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count only carbon atoms
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == carbonyl_c_idx:
                    continue
                if neighbor_idx not in visited_acyl and neighbor_idx not in CoA_atom_indices:
                    traverse_acyl_chain(neighbor_idx)
        traverse_acyl_chain(alpha_c_idx)

        # Check if carbon count is within 13 to 22
        if 13 <= carbon_count <= 22:
            return True, f"Contains long-chain fatty acyl group with {carbon_count} carbons"
        else:
            continue  # Try next thioester linkage if any

    # If no suitable acyl chain found
    return False, "No long-chain fatty acyl chain of length 13-22 carbons found"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57395',
        'name': 'long-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.',
        'parents': ['CHEBI:37554']
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
    'stdout': None
}