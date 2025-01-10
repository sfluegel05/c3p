"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57395 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the condensation of the thiol group of coenzyme A with
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

    # Define the Coenzyme A (CoA) molecule from SMILES
    coa_smiles = "C[C@H](O)[C@H](OP(=O)(O)OCC1OC(O)[C@H](OP(=O)(O)O)[C@@H](O)[C@H]1O)n1cnc2c(N)ncnc12"
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Failed to construct CoA molecule"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"

    # Define thioester linkage pattern (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if thioester_pattern is None:
        return False, "Failed to construct thioester pattern"

    # Find the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify the fatty acyl chain attached via thioester linkage
    for match in thioester_matches:
        sulfur_idx = match[0]
        carbonyl_c_idx = match[1]

        # Get the carbon adjacent to the carbonyl carbon (start of fatty acyl chain)
        fatty_acyl_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        neighbors = [atom.GetIdx() for atom in fatty_acyl_atom.GetNeighbors() if atom.GetIdx() != sulfur_idx]
        if not neighbors:
            continue
        fatty_acyl_start_idx = neighbors[0]

        # Use a BFS traversal to find the length of the fatty acyl chain
        visited = set()
        queue = [fatty_acyl_start_idx]
        carbon_count = 0

        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                carbon_count += 1
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                        queue.append(neighbor_idx)

        # Check if carbon count is within the long-chain fatty acid range (C13 to C22)
        if 13 <= carbon_count <= 22:
            return True, f"Contains long-chain fatty acyl group with {carbon_count} carbons"
        else:
            return False, f"Fatty acyl chain length is {carbon_count} carbons, not in range 13-22"

    return False, "Failed to identify long-chain fatty acyl chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57395',
        'name': 'long-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:64479']
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
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}